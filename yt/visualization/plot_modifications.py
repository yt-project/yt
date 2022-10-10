import inspect
import re
import warnings
from abc import ABC, abstractmethod
from functools import update_wrapper
from numbers import Integral, Number, Real
from typing import Any, Dict, Optional, Tuple, Type, Union

import matplotlib
import numpy as np

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.data_objects.data_containers import YTDataContainer
from yt.data_objects.level_sets.clump_handling import Clump
from yt.data_objects.selection_objects.cut_region import YTCutRegion
from yt.data_objects.static_output import Dataset
from yt.frontends.ytdata.data_structures import YTClumpContainer
from yt.funcs import is_sequence, mylog, validate_width_tuple
from yt.geometry.geometry_handler import is_curvilinear
from yt.geometry.unstructured_mesh_handler import UnstructuredIndex
from yt.units import dimensions
from yt.units.yt_array import YTArray, YTQuantity, uhstack  # type: ignore
from yt.utilities.exceptions import YTDataTypeUnsupported, YTUnsupportedPlotCallback
from yt.utilities.lib.geometry_utils import triangle_plane_intersect
from yt.utilities.lib.line_integral_convolution import line_integral_convolution_2d
from yt.utilities.lib.mesh_triangulation import triangulate_indices
from yt.utilities.lib.pixelization_routines import (
    pixelize_cartesian,
    pixelize_off_axis_cartesian,
)
from yt.utilities.math_utils import periodic_ray
from yt.utilities.on_demand_imports import NotAModule
from yt.visualization._commons import (
    _swap_arg_pair_order,
    _swap_axes_extents,
    invalidate_plot,
)
from yt.visualization.base_plot_types import CallbackWrapper
from yt.visualization.image_writer import apply_colormap
from yt.visualization.plot_window import PWViewerMPL

callback_registry: Dict[str, Type["PlotCallback"]] = {}


def _validate_factor_tuple(factor) -> Tuple[int, int]:
    if (
        is_sequence(factor)
        and len(factor) == 2
        and all(isinstance(_, Integral) for _ in factor)
    ):
        # - checking for "is_sequence" allows lists, numpy arrays and other containers
        # - checking for Integral type allows numpy integer types
        # in any case we return a with strict typing
        return (int(factor[0]), int(factor[1]))
    elif isinstance(factor, Integral):
        return (int(factor), int(factor))
    else:
        raise TypeError(
            f"Expected a single, or a pair of integers, received {factor!r}"
        )


class PlotCallback(ABC):
    # _supported_geometries is set by subclasses of PlotCallback to a tuple of
    # strings corresponding to the names of the geometries that a callback
    # supports.  By default it is None, which means it supports everything.
    # Note that if there's a coord_system parameter that is set to "axis" or
    # "figure" this is disregarded.  If "force" is included in the tuple, it
    # will *not* check whether or not the coord_system is in axis or figure,
    # and will only look at the geometries.
    _supported_geometries: Optional[Tuple[str, ...]] = None
    _incompatible_plot_types: Tuple[str, ...] = tuple()

    def __init_subclass__(cls, *args, **kwargs):
        if inspect.isabstract(cls):
            return

        # register class
        callback_registry[cls.__name__] = cls

        # create a PWViewerMPL method by wrapping __init__
        if cls.__init__.__doc__ is None:
            # allow docstring definition at the class level instead of __init__
            cls.__init__.__doc__ = cls.__doc__

        supported_geometries = cls._supported_geometries
        incompatible_plot_types = cls._incompatible_plot_types
        type_name = cls._type_name

        @invalidate_plot
        def closure(self, *args, **kwargs):
            nonlocal supported_geometries
            nonlocal incompatible_plot_types
            nonlocal type_name

            geom = self.ds.geometry
            if not (
                supported_geometries is None
                or geom in supported_geometries
                or (
                    kwargs.get("coord_system") in ("axis", "figure")
                    and "force" not in supported_geometries
                )
            ):
                raise YTDataTypeUnsupported(geom, supported_geometries)
            if self._plot_type in incompatible_plot_types:
                raise YTUnsupportedPlotCallback(type_name, self._plot_type)
            self._callbacks.append(cls(*args, **kwargs))
            return self

        update_wrapper(
            wrapper=closure,
            wrapped=cls.__init__,
            assigned=("__annotations__", "__doc__"),
        )

        method_name = "annotate_" + type_name
        closure.__name__ = method_name
        setattr(PWViewerMPL, method_name, closure)

    @abstractmethod
    def __init__(self, *args, **kwargs) -> None:
        pass

    @abstractmethod
    def __call__(self, plot: CallbackWrapper) -> None:
        pass

    def _project_coords(self, plot, coord):
        """
        Convert coordinates from simulation data coordinates to projected
        data coordinates.  Simulation data coordinates are three dimensional,
        and can either be specified as a YTArray or as a list or array in
        code_length units.  Projected data units are 2D versions of the
        simulation data units relative to the axes of the final plot.
        """
        if len(coord) == 3:
            if not isinstance(coord, YTArray):
                coord_copy = plot.data.ds.arr(coord, "code_length")
            else:
                # coord is being copied so that if the user has a unyt_array already
                # we don't change the user's version
                coord_copy = coord.to("code_length")
            ax = plot.data.axis
            # if this is an on-axis projection or slice, then
            # just grab the appropriate 2 coords for the on-axis view
            if ax >= 0 and ax <= 2:
                (xi, yi) = (
                    plot.data.ds.coordinates.x_axis[ax],
                    plot.data.ds.coordinates.y_axis[ax],
                )
                ret_coord = (coord_copy[xi], coord_copy[yi])

            # if this is an off-axis project or slice (ie cutting plane)
            # we have to calculate where the data coords fall in the projected
            # plane
            elif ax == 4:
                # transpose is just to get [[x1,x2,...],[y1,y2,...],[z1,z2,...]]
                # in the same order as plot.data.center for array arithmetic
                coord_vectors = coord_copy.transpose() - plot.data.center
                x = np.dot(coord_vectors, plot.data.orienter.unit_vectors[1])
                y = np.dot(coord_vectors, plot.data.orienter.unit_vectors[0])
                # Transpose into image coords. Due to VR being not a
                # right-handed coord system
                ret_coord = (y, x)
            else:
                raise ValueError("Object being plotted must have a `data.axis` defined")

        # if the position is already two-coords, it is expected to be
        # in the proper projected orientation
        else:
            raise ValueError("'data' coordinates must be 3 dimensions")
        return ret_coord

    def _convert_to_plot(self, plot, coord, offset=True):
        """
        Convert coordinates from projected data coordinates to PlotWindow
        plot coordinates.  Projected data coordinates are two dimensional
        and refer to the location relative to the specific axes being plotted,
        although still in simulation units.  PlotWindow plot coordinates
        are locations as found in the final plot, usually with the origin
        in the center of the image and the extent of the image defined by
        the final plot axis markers.
        """
        # coord should be a 2 x ncoord array-like datatype.
        try:
            ncoord = np.array(coord).shape[1]
        except IndexError:
            ncoord = 1

        # Convert the data and plot limits to tiled numpy arrays so that
        # convert_to_plot is automatically vectorized.
        phy_bounds = self._physical_bounds(plot)
        x0, x1, y0, y1 = (np.array(np.tile(xyi, ncoord)) for xyi in phy_bounds)
        plt_bounds = self._plot_bounds(plot)
        xx0, xx1, yy0, yy1 = (np.tile(xyi, ncoord) for xyi in plt_bounds)

        try:
            ccoord = np.array(coord.to("code_length"))
        except AttributeError:
            ccoord = np.array(coord)

        # We need a special case for when we are only given one coordinate.
        if ccoord.shape == (2,):
            return np.array(
                [
                    ((ccoord[0] - x0) / (x1 - x0) * (xx1 - xx0) + xx0)[0],
                    ((ccoord[1] - y0) / (y1 - y0) * (yy1 - yy0) + yy0)[0],
                ]
            )
        else:
            return np.array(
                [
                    (ccoord[0][:] - x0) / (x1 - x0) * (xx1 - xx0) + xx0,
                    (ccoord[1][:] - y0) / (y1 - y0) * (yy1 - yy0) + yy0,
                ]
            )

    def _sanitize_coord_system(self, plot, coord, coord_system):
        """
        Given a set of one or more x,y (and z) coordinates and a coordinate
        system, convert the coordinates (and transformation) ready for final
        plotting.

        Parameters
        ----------

        plot: a PlotMPL subclass
           The plot that we are converting coordinates for

        coord: array-like
           Coordinates in some coordinate system: [x,y,z].
           Alternatively, can specify multiple coordinates as:
           [[x1,x2,...,xn], [y1, y2,...,yn], [z1,z2,...,zn]]

        coord_system: string

            Possible values include:

            * ``'data'``
                3D data coordinates relative to original dataset

            * ``'plot'``
                2D coordinates as defined by the final axis locations

            * ``'axis'``
                2D coordinates within the axis object from (0,0) in lower left
                to (1,1) in upper right.  Same as matplotlib axis coords.

            * ``'figure'``
                2D coordinates within figure object from (0,0) in lower left
                to (1,1) in upper right.  Same as matplotlib figure coords.
        """
        # Assure coords are either a YTArray or numpy array
        coord = np.asanyarray(coord, dtype="float64")
        # if in data coords, project them to plot coords
        if coord_system == "data":
            if len(coord) < 3:
                raise ValueError(
                    "Coordinates in 'data' coordinate system need to be in 3D"
                )
            coord = self._project_coords(plot, coord)
            coord = self._convert_to_plot(plot, coord)
        # if in plot coords, define the transform correctly
        if coord_system == "data" or coord_system == "plot":
            self.transform = plot._axes.transData
            return coord
        # if in axis coords, define the transform correctly
        if coord_system == "axis":
            self.transform = plot._axes.transAxes
            if len(coord) > 2:
                raise ValueError(
                    "Coordinates in 'axis' coordinate system need to be in 2D"
                )
            return coord
        # if in figure coords, define the transform correctly
        elif coord_system == "figure":
            self.transform = plot._figure.transFigure
            return coord
        else:
            raise ValueError(
                "Argument coord_system must have a value of "
                "'data', 'plot', 'axis', or 'figure'."
            )

    def _physical_bounds(self, plot):
        xlims = tuple(v.in_units("code_length") for v in plot.xlim)
        ylims = tuple(v.in_units("code_length") for v in plot.ylim)
        # _swap_axes note: do NOT need to unswap here because plot (a CallbackWrapper
        # instance) stores the x, y lims of the underlying data object.
        return xlims + ylims

    def _plot_bounds(self, plot):
        xlims = plot._axes.get_xlim()
        ylims = plot._axes.get_ylim()
        # _swap_axes note: because we are getting the plot limits from the axes
        # object, if the axes have been swapped, these will be reversed from the
        # _physical_bounds. So we need to unswap here, but not in _physical_bounds.
        if plot._swap_axes:
            return ylims + xlims
        return xlims + ylims

    def _pixel_scale(self, plot):
        x0, x1, y0, y1 = self._physical_bounds(plot)
        xx0, xx1, yy0, yy1 = self._plot_bounds(plot)
        dx = (xx1 - xx0) / (x1 - x0)
        dy = (yy1 - yy0) / (y1 - y0)
        return dx, dy

    def _set_font_properties(self, plot, labels, **kwargs):
        """
        This sets all of the text instances created by a callback to have
        the same font size and properties as all of the other fonts in the
        figure.  If kwargs are set, they override the defaults.
        """
        # This is a little messy because there is no trivial way to update
        # a MPL.font_manager.FontProperties object with new attributes
        # aside from setting them individually.  So we pick out the relevant
        # MPL.Text() kwargs from the local kwargs and let them override the
        # defaults.
        local_font_properties = plot.font_properties.copy()

        # Turn off the default TT font file, otherwise none of this works.
        local_font_properties.set_file(None)
        local_font_properties.set_family("stixgeneral")

        if "family" in kwargs:
            local_font_properties.set_family(kwargs["family"])
        if "file" in kwargs:
            local_font_properties.set_file(kwargs["file"])
        if "fontconfig_pattern" in kwargs:
            local_font_properties.set_fontconfig_pattern(kwargs["fontconfig_pattern"])
        if "name" in kwargs:
            local_font_properties.set_name(kwargs["name"])
        if "size" in kwargs:
            local_font_properties.set_size(kwargs["size"])
        if "slant" in kwargs:
            local_font_properties.set_slant(kwargs["slant"])
        if "stretch" in kwargs:
            local_font_properties.set_stretch(kwargs["stretch"])
        if "style" in kwargs:
            local_font_properties.set_style(kwargs["style"])
        if "variant" in kwargs:
            local_font_properties.set_variant(kwargs["variant"])
        if "weight" in kwargs:
            local_font_properties.set_weight(kwargs["weight"])

        # For each label, set the font properties and color to the figure
        # defaults if not already set in the callback itself
        for label in labels:
            if plot.font_color is not None and "color" not in kwargs:
                label.set_color(plot.font_color)
            label.set_fontproperties(local_font_properties)

    def _set_plot_limits(self, plot, extent=None) -> None:
        """
        calls set_xlim, set_ylim for plot, accounting for swapped axes

        Parameters
        ----------
        plot : CallbackWrapper
            a CallbackWrapper instance
        extent : tuple or list
            The raw extent (prior to swapping). if None, will fetch it.
        """
        if extent is None:
            extent = self._plot_bounds(plot)

        if plot._swap_axes:
            extent = _swap_axes_extents(extent)

        plot._axes.set_xlim(extent[0], extent[1])
        plot._axes.set_ylim(extent[2], extent[3])

    @staticmethod
    def _sanitize_xy_order(plot, *args):
        """
        flips x-y pairs of plot arguments if needed

        Parameters
        ----------
        plot : CallbackWrapper
            a CallbackWrapper instance
        *args
            x, y plot arguments, must have an even number of *args

        Returns
        -------
        tuple
            either the original args or new args, with (x, y) pairs switched. i.e.,

            _sanitize_xy_order(plot, x, y, px, py) returns:
                x, y, px, py if plot._swap_axes is False
                y, x, py, px if plot._swap_axes is True

        """
        if plot._swap_axes:
            return _swap_arg_pair_order(*args)
        return args


class VelocityCallback(PlotCallback):
    """
    Adds a 'quiver' plot of velocity to the plot, skipping all but
    every *factor* datapoint. *scale* is the data units per arrow
    length unit using *scale_units* and *plot_args* allows you to
    pass in matplotlib arguments (see matplotlib.axes.Axes.quiver
    for more info). if *normalize* is True, the velocity fields
    will be scaled by their local (in-plane) length, allowing
    morphological features to be more clearly seen for fields
    with substantial variation in field strength.
    """

    _type_name = "velocity"
    _supported_geometries = (
        "cartesian",
        "spectral_cube",
        "polar",
        "cylindrical",
        "spherical",
    )
    _incompatible_plot_types = ("OffAxisProjection", "Particle")

    def __init__(
        self,
        factor: Union[Tuple[int, int], int] = 16,
        *,
        scale=None,
        scale_units=None,
        normalize=False,
        plot_args=None,
        **kwargs,
    ):
        self.factor = _validate_factor_tuple(factor)
        self.scale = scale
        self.scale_units = scale_units
        self.normalize = normalize
        if plot_args is not None:
            issue_deprecation_warning(
                "`plot_args` is deprecated. "
                "You can now pass arbitrary keyword arguments instead of a dictionary.",
                since="4.1.0",
            )
            plot_args.update(kwargs)
        else:
            plot_args = kwargs

        self.plot_args = plot_args

    def __call__(self, plot):
        ftype = plot.data._current_fluid_type
        # Instantiation of these is cheap
        if plot._type_name == "CuttingPlane":
            if is_curvilinear(plot.data.ds.geometry):
                raise NotImplementedError(
                    "Velocity annotation for cutting \
                    plane is not supported for %s geometry"
                    % plot.data.ds.geometry
                )
        if plot._type_name == "CuttingPlane":
            qcb = CuttingQuiverCallback(
                (ftype, "cutting_plane_velocity_x"),
                (ftype, "cutting_plane_velocity_y"),
                factor=self.factor,
                scale=self.scale,
                normalize=self.normalize,
                scale_units=self.scale_units,
                **self.plot_args,
            )
        else:
            xax = plot.data.ds.coordinates.x_axis[plot.data.axis]
            yax = plot.data.ds.coordinates.y_axis[plot.data.axis]
            axis_names = plot.data.ds.coordinates.axis_name

            bv = plot.data.get_field_parameter("bulk_velocity")
            if bv is not None:
                bv_x = bv[xax]
                bv_y = bv[yax]
            else:
                bv_x = bv_y = 0

            if (
                plot.data.ds.geometry in ["polar", "cylindrical"]
                and axis_names[plot.data.axis] == "z"
            ):
                # polar_z and cyl_z is aligned with carteian_z
                # should convert r-theta plane to x-y plane
                xv = (ftype, "velocity_cartesian_x")
                yv = (ftype, "velocity_cartesian_y")
            elif plot.data.ds.geometry == "spherical":
                if axis_names[plot.data.axis] == "phi":
                    xv = (ftype, "velocity_cylindrical_radius")
                    yv = (ftype, "velocity_cylindrical_z")
                elif axis_names[plot.data.axis] == "theta":
                    xv = (ftype, "velocity_conic_x")
                    yv = (ftype, "velocity_conic_y")
                else:
                    raise NotImplementedError(
                        f"annotate_velocity is missing support for normal={axis_names[plot.data.axis]!r}"
                    )
            else:
                # for other cases (even for cylindrical geometry),
                # orthogonal planes are generically Cartesian
                xv = (ftype, f"velocity_{axis_names[xax]}")
                yv = (ftype, f"velocity_{axis_names[yax]}")

            # determine the full fields including field type
            xv = plot.data._determine_fields(xv)[0]
            yv = plot.data._determine_fields(yv)[0]

            qcb = QuiverCallback(
                xv,
                yv,
                factor=self.factor,
                scale=self.scale,
                scale_units=self.scale_units,
                normalize=self.normalize,
                bv_x=bv_x,
                bv_y=bv_y,
                **self.plot_args,
            )
        return qcb(plot)


class MagFieldCallback(PlotCallback):
    """
    Adds a 'quiver' plot of magnetic field to the plot, skipping all but
    every *factor* datapoint. *scale* is the data units per arrow
    length unit using *scale_units* and *plot_args* allows you to pass
    in matplotlib arguments (see matplotlib.axes.Axes.quiver for more info).
    if *normalize* is True, the magnetic fields will be scaled by their
    local (in-plane) length, allowing morphological features to be more
    clearly seen for fields with substantial variation in field strength.
    """

    _type_name = "magnetic_field"
    _supported_geometries = (
        "cartesian",
        "spectral_cube",
        "polar",
        "cylindrical",
        "spherical",
    )
    _incompatible_plot_types = ("OffAxisProjection", "Particle")

    def __init__(
        self,
        factor: Union[Tuple[int, int], int] = 16,
        *,
        scale=None,
        scale_units=None,
        normalize=False,
        plot_args=None,
        **kwargs,
    ):
        self.factor = _validate_factor_tuple(factor)
        self.scale = scale
        self.scale_units = scale_units
        self.normalize = normalize
        if plot_args is not None:
            issue_deprecation_warning(
                "`plot_args` is deprecated. "
                "You can now pass arbitrary keyword arguments instead of a dictionary.",
                since="4.1.0",
            )
            plot_args.update(kwargs)
        else:
            plot_args = kwargs
        self.plot_args = plot_args

    def __call__(self, plot):
        ftype = plot.data._current_fluid_type
        # Instantiation of these is cheap
        if plot._type_name == "CuttingPlane":
            if is_curvilinear(plot.data.ds.geometry):
                raise NotImplementedError(
                    "Magnetic field annotation for cutting \
                    plane is not supported for %s geometry"
                    % plot.data.ds.geometry
                )
            qcb = CuttingQuiverCallback(
                (ftype, "cutting_plane_magnetic_field_x"),
                (ftype, "cutting_plane_magnetic_field_y"),
                factor=self.factor,
                scale=self.scale,
                scale_units=self.scale_units,
                normalize=self.normalize,
                **self.plot_args,
            )
        else:
            xax = plot.data.ds.coordinates.x_axis[plot.data.axis]
            yax = plot.data.ds.coordinates.y_axis[plot.data.axis]
            axis_names = plot.data.ds.coordinates.axis_name

            if (
                plot.data.ds.geometry in ["polar", "cylindrical"]
                and axis_names[plot.data.axis] == "z"
            ):
                # polar_z and cyl_z is aligned with carteian_z
                # should convert r-theta plane to x-y plane
                xv = (ftype, "magnetic_field_cartesian_x")
                yv = (ftype, "magnetic_field_cartesian_y")
            elif plot.data.ds.geometry == "spherical":
                if axis_names[plot.data.axis] == "phi":
                    xv = (ftype, "magnetic_field_cylindrical_radius")
                    yv = (ftype, "magnetic_field_cylindrical_z")
                elif axis_names[plot.data.axis] == "theta":
                    xv = (ftype, "magnetic_field_conic_x")
                    yv = (ftype, "magnetic_field_conic_y")
                else:
                    raise NotImplementedError(
                        f"annotate_magnetic_field is missing support for normal={axis_names[plot.data.axis]!r}"
                    )
            else:
                # for other cases (even for cylindrical geometry),
                # orthogonal planes are generically Cartesian
                xv = (ftype, f"magnetic_field_{axis_names[xax]}")
                yv = (ftype, f"magnetic_field_{axis_names[yax]}")

            qcb = QuiverCallback(
                xv,
                yv,
                factor=self.factor,
                scale=self.scale,
                scale_units=self.scale_units,
                normalize=self.normalize,
                **self.plot_args,
            )
        return qcb(plot)


class BaseQuiverCallback(PlotCallback, ABC):
    _incompatible_plot_types = ("OffAxisProjection", "Particle")

    def __init__(
        self,
        field_x,
        field_y,
        field_c=None,
        *,
        factor: Union[Tuple[int, int], int] = 16,
        scale=None,
        scale_units=None,
        normalize=False,
        plot_args=None,
        **kwargs,
    ):
        self.field_x = field_x
        self.field_y = field_y
        self.field_c = field_c
        self.factor = _validate_factor_tuple(factor)
        self.scale = scale
        self.scale_units = scale_units
        self.normalize = normalize
        if plot_args is None:
            plot_args = kwargs
        else:
            issue_deprecation_warning(
                "`plot_args` is deprecated. "
                "You can now pass arbitrary keyword arguments instead of a dictionary.",
                since="4.1.0",
            )
            plot_args.update(kwargs)

        self.plot_args = plot_args

    @abstractmethod
    def _get_quiver_data(self, plot, bounds: tuple, nx: int, ny: int):
        # must return (pixX, pixY, pixC) arrays (pixC can be None)
        pass

    def __call__(self, plot):

        # construct mesh
        bounds = self._physical_bounds(plot)
        nx = plot.raw_image_shape[1] // self.factor[0]
        ny = plot.raw_image_shape[0] // self.factor[1]
        xx0, xx1, yy0, yy1 = self._plot_bounds(plot)
        X, Y = np.meshgrid(
            np.linspace(xx0, xx1, nx, endpoint=True),
            np.linspace(yy0, yy1, ny, endpoint=True),
        )

        pixX, pixY, pixC = self._get_quiver_data(plot, bounds, nx, ny)

        retv = self._finalize(plot, X, Y, pixX, pixY, pixC)
        self._set_plot_limits(plot, (xx0, xx1, yy0, yy1))
        return retv

    def _finalize(self, plot, X, Y, pixX, pixY, pixC):
        if self.normalize:
            nn = np.sqrt(pixX**2 + pixY**2)
            nn = np.where(nn == 0, 1, nn)
            pixX /= nn
            pixY /= nn

        X, Y, pixX, pixY = self._sanitize_xy_order(plot, X, Y, pixX, pixY)
        # quiver plots ignore x, y axes inversions when using angles="uv" (the
        # default), so reverse the direction of the vectors when flipping the axis.
        if self.plot_args.get("angles", None) != "xy":
            if plot._flip_vertical:
                pixY = -1 * pixY
            if plot._flip_horizontal:
                pixX = -1 * pixX

        args = [X, Y, pixX, pixY]
        if pixC is not None:
            args.append(pixC)

        kwargs = dict(
            scale=self.scale,
            scale_units=self.scale_units,
        )
        kwargs.update(self.plot_args)
        return plot._axes.quiver(*args, **kwargs)


class QuiverCallback(BaseQuiverCallback):
    """
    Adds a 'quiver' plot to any plot, using the *field_x* and *field_y*
    from the associated data, skipping every *factor* pixels.
    *field_c* is an optional field name used for color.
    *scale* is the data units per arrow length unit using *scale_units*
    and *plot_args* allows you to pass in matplotlib arguments (see
    matplotlib.axes.Axes.quiver for more info). if *normalize* is True,
    the fields will be scaled by their local (in-plane) length, allowing
    morphological features to be more clearly seen for fields with
    substantial variation in field strength.
    """

    _type_name = "quiver"
    _supported_geometries: Tuple[str, ...] = (
        "cartesian",
        "spectral_cube",
        "polar",
        "cylindrical",
        "spherical",
    )

    def __init__(
        self,
        field_x,
        field_y,
        field_c=None,
        *,
        factor: Union[Tuple[int, int], int] = 16,
        scale=None,
        scale_units=None,
        normalize=False,
        bv_x=0,
        bv_y=0,
        plot_args=None,
        **kwargs,
    ):
        super().__init__(
            field_x,
            field_y,
            field_c,
            factor=factor,
            scale=scale,
            scale_units=scale_units,
            normalize=normalize,
            plot_args=plot_args,
            **kwargs,
        )
        self.bv_x = bv_x
        self.bv_y = bv_y

    def _get_quiver_data(self, plot, bounds: tuple, nx: int, ny: int):
        # calls the pixelizer, returns pixX, pixY, pixC arrays

        def transform(field_name, vector_value):
            field_units = plot.data[field_name].units

            def _transformed_field(field, data):
                return data[field_name] - data.ds.arr(vector_value, field_units)

            plot.data.ds.add_field(
                ("gas", f"transformed_{field_name}"),
                sampling_type="cell",
                function=_transformed_field,
                units=field_units,
                display_field=False,
            )

        if self.bv_x != 0.0 or self.bv_x != 0.0:
            # We create a relative vector field
            transform(self.field_x, self.bv_x)
            transform(self.field_y, self.bv_y)
            field_x = f"transformed_{self.field_x}"
            field_y = f"transformed_{self.field_y}"
        else:
            field_x, field_y = self.field_x, self.field_y

        periodic = int(any(plot.data.ds.periodicity))
        pixX = plot.data.ds.coordinates.pixelize(
            plot.data.axis,
            plot.data,
            field_x,
            bounds,
            (nx, ny),
            False,  # antialias
            periodic,
        )
        pixY = plot.data.ds.coordinates.pixelize(
            plot.data.axis,
            plot.data,
            field_y,
            bounds,
            (nx, ny),
            False,  # antialias
            periodic,
        )

        if self.field_c is not None:
            pixC = plot.data.ds.coordinates.pixelize(
                plot.data.axis,
                plot.data,
                self.field_c,
                bounds,
                (nx, ny),
                False,  # antialias
                periodic,
            )
        else:
            pixC = None

        return pixX, pixY, pixC


class ContourCallback(PlotCallback):
    """
    Add contours in *field* to the plot. *levels* governs the number of
    contours generated, *factor* governs the number of points used in the
    interpolation, *take_log* governs how it is contoured and *clim* gives
    the (upper, lower) limits for contouring.  An alternate data source can be
    specified with *data_source*, but by default the plot's data source will be
    queried.
    """

    _type_name = "contour"
    _supported_geometries = ("cartesian", "spectral_cube", "cylindrical")
    _incompatible_plot_types = ("Particle",)

    def __init__(
        self,
        field: Union[Tuple[str, str], str],
        levels: int = 5,
        *,
        factor: Union[Tuple[int, int], int] = 4,
        clim: Optional[Tuple[float, float]] = None,
        label: bool = False,
        take_log: Optional[bool] = None,
        data_source: YTDataContainer = None,
        plot_args: Optional[Dict[str, Any]] = None,
        label_args: Optional[Dict[str, Any]] = None,
        text_args: Optional[Dict[str, Any]] = None,
        ncont: Optional[int] = None,  # deprecated
    ) -> None:
        def_plot_args = {"colors": "k", "linestyles": "solid"}
        def_text_args = {"colors": "w"}
        if ncont is not None:
            issue_deprecation_warning(
                "The `ncont` keyword argument is deprecated, use `levels` instead.",
                since="4.1.0",
            )
            levels = ncont
        if clim is not None and not isinstance(levels, (int, np.integer)):
            raise TypeError(f"clim requires levels be an integer, received {levels}")

        self.levels = levels
        self.field = field
        self.factor = _validate_factor_tuple(factor)
        self.clim = clim
        self.take_log = take_log
        if plot_args is None:
            plot_args = def_plot_args
        self.plot_args = plot_args
        self.label = label
        if label_args is not None:
            text_args = label_args
            issue_deprecation_warning(
                "The label_args keyword is deprecated.  Please use "
                "the text_args keyword instead.",
                since="3.2",
                removal="4.2",
            )
        if text_args is None:
            text_args = def_text_args
        self.text_args = text_args
        self.data_source = data_source

    def __call__(self, plot):
        from matplotlib.tri import LinearTriInterpolator, Triangulation

        # These need to be in code_length
        x0, x1, y0, y1 = self._physical_bounds(plot)

        # These are in plot coordinates, which may not be code coordinates.
        xx0, xx1, yy0, yy1 = self._plot_bounds(plot)

        # See the note about rows/columns in the pixelizer for more information
        # on why we choose the bounds we do
        numPoints_x = plot.raw_image_shape[1]
        numPoints_y = plot.raw_image_shape[0]

        # Go from data->plot coordinates
        dx, dy = self._pixel_scale(plot)

        # We want xi, yi in plot coordinates
        xi, yi = np.mgrid[
            xx0 : xx1 : numPoints_x / (self.factor[0] * 1j),
            yy0 : yy1 : numPoints_y / (self.factor[1] * 1j),
        ]
        data = self.data_source or plot.data

        if plot._type_name in ["CuttingPlane", "Projection", "Slice"]:
            if plot._type_name == "CuttingPlane":
                x = data["px"] * dx
                y = data["py"] * dy
                z = data[self.field]
            elif plot._type_name in ["Projection", "Slice"]:
                # Makes a copy of the position fields "px" and "py" and adds the
                # appropriate shift to the copied field.

                AllX = np.zeros(data["px"].size, dtype="bool")
                AllY = np.zeros(data["py"].size, dtype="bool")
                XShifted = data["px"].copy()
                YShifted = data["py"].copy()
                dom_x, dom_y = plot._period
                for shift in np.mgrid[-1:1:3j]:
                    xlim = (data["px"] + shift * dom_x >= x0) & (
                        data["px"] + shift * dom_x <= x1
                    )
                    ylim = (data["py"] + shift * dom_y >= y0) & (
                        data["py"] + shift * dom_y <= y1
                    )
                    XShifted[xlim] += shift * dom_x
                    YShifted[ylim] += shift * dom_y
                    AllX |= xlim
                    AllY |= ylim

                # At this point XShifted and YShifted are the shifted arrays of
                # position data in data coordinates
                wI = AllX & AllY

                # This converts XShifted and YShifted into plot coordinates
                x = ((XShifted[wI] - x0) * dx).ndarray_view() + xx0
                y = ((YShifted[wI] - y0) * dy).ndarray_view() + yy0
                z = data[self.field][wI]

            # Both the input and output from the triangulator are in plot
            # coordinates
            triangulation = Triangulation(x, y)
            zi = LinearTriInterpolator(triangulation, z)(xi, yi)
        elif plot._type_name == "OffAxisProjection":
            zi = plot.frb[self.field][:: self.factor[0], :: self.factor[1]].transpose()

        take_log: bool
        if self.take_log is not None:
            take_log = self.take_log
        else:
            field = data._determine_fields([self.field])[0]
            take_log = plot.ds._get_field_info(*field).take_log

        if take_log:
            zi = np.log10(zi)

        clim: Union[Tuple[Real, Real], None]
        if take_log and self.clim is not None:
            clim = np.log10(self.clim[0]), np.log10(self.clim[1])
        else:
            clim = self.clim

        levels: Union[np.ndarray, int]
        if clim is not None:
            levels = np.linspace(clim[0], clim[1], self.levels)
        else:
            levels = self.levels

        xi, yi = self._sanitize_xy_order(plot, xi, yi)
        cset = plot._axes.contour(xi, yi, zi, levels, **self.plot_args)
        self._set_plot_limits(plot, (xx0, xx1, yy0, yy1))

        if self.label:
            plot._axes.clabel(cset, **self.text_args)


class GridBoundaryCallback(PlotCallback):
    """
    Draws grids on an existing PlotWindow object.  Adds grid boundaries to a
    plot, optionally with alpha-blending. By default, colors different levels of
    grids with different colors going from white to black, but you can change to
    any arbitrary colormap with cmap keyword, to all black grid edges for all
    levels with cmap=None and edgecolors=None, or to an arbitrary single color
    for grid edges with edgecolors='YourChosenColor' defined in any of the
    standard ways (e.g., edgecolors='white', edgecolors='r',
    edgecolors='#00FFFF', or edgecolor='0.3', where the last is a float in 0-1
    scale indicating gray).  Note that setting edgecolors overrides cmap if you
    have both set to non-None values.  Cutoff for display is at min_pix
    wide. draw_ids puts the grid id a the corner of the grid (but its not so
    great in projections...).  id_loc determines which corner holds the grid id.
    One can set min and maximum level of grids to display, and
    can change the linewidth of the displayed grids.
    """

    _type_name = "grids"
    _supported_geometries = ("cartesian", "spectral_cube", "cylindrical")
    _incompatible_plot_types = ("OffAxisSlice", "OffAxisProjection", "Particle")

    def __init__(
        self,
        alpha=0.7,
        min_pix=1,
        min_pix_ids=20,
        draw_ids=False,
        id_loc=None,
        periodic=True,
        min_level=None,
        max_level=None,
        cmap="B-W LINEAR_r",
        edgecolors=None,
        linewidth=1.0,
    ):
        self.alpha = alpha
        self.min_pix = min_pix
        self.min_pix_ids = min_pix_ids
        self.draw_ids = draw_ids  # put grid numbers in the corner.
        if id_loc is None:
            self.id_loc = "lower left"
        else:
            self.id_loc = id_loc.lower()  # Make case-insensitive
            if not self.draw_ids:
                mylog.warning(
                    "Supplied id_loc but draw_ids is False. Not drawing grid ids"
                )
        self.periodic = periodic
        self.min_level = min_level
        self.max_level = max_level
        self.linewidth = linewidth
        self.cmap = cmap
        self.edgecolors = edgecolors

    def __call__(self, plot):
        if plot.data.ds.geometry == "cylindrical" and plot.data.ds.dimensionality == 3:
            raise NotImplementedError(
                "Grid annotation is only supported for \
                for 2D cylindrical geometry, not 3D"
            )
        from matplotlib.colors import colorConverter

        x0, x1, y0, y1 = self._physical_bounds(plot)
        xx0, xx1, yy0, yy1 = self._plot_bounds(plot)
        (dx, dy) = self._pixel_scale(plot)
        (ypix, xpix) = plot.raw_image_shape
        ax = plot.data.axis
        px_index = plot.data.ds.coordinates.x_axis[ax]
        py_index = plot.data.ds.coordinates.y_axis[ax]
        DW = plot.data.ds.domain_width
        if self.periodic:
            pxs, pys = np.mgrid[-1:1:3j, -1:1:3j]
        else:
            pxs, pys = np.mgrid[0:0:1j, 0:0:1j]
        GLE, GRE, levels, block_ids = [], [], [], []
        for block, _mask in plot.data.blocks:
            GLE.append(block.LeftEdge.in_units("code_length"))
            GRE.append(block.RightEdge.in_units("code_length"))
            levels.append(block.Level)
            block_ids.append(block.id)
        if len(GLE) == 0:
            return
        # Retain both units and registry
        GLE = plot.ds.arr(GLE, units=GLE[0].units)
        GRE = plot.ds.arr(GRE, units=GRE[0].units)
        levels = np.array(levels)
        min_level = self.min_level or 0
        max_level = self.max_level or levels.max()

        # sort the four arrays in order of ascending level, this makes images look nicer
        new_indices = np.argsort(levels)
        levels = levels[new_indices]
        GLE = GLE[new_indices]
        GRE = GRE[new_indices]
        block_ids = np.array(block_ids)[new_indices]

        for px_off, py_off in zip(pxs.ravel(), pys.ravel()):
            pxo = px_off * DW[px_index]
            pyo = py_off * DW[py_index]
            left_edge_x = np.array((GLE[:, px_index] + pxo - x0) * dx) + xx0
            left_edge_y = np.array((GLE[:, py_index] + pyo - y0) * dy) + yy0
            right_edge_x = np.array((GRE[:, px_index] + pxo - x0) * dx) + xx0
            right_edge_y = np.array((GRE[:, py_index] + pyo - y0) * dy) + yy0
            xwidth = xpix * (right_edge_x - left_edge_x) / (xx1 - xx0)
            ywidth = ypix * (right_edge_y - left_edge_y) / (yy1 - yy0)
            visible = np.logical_and(
                np.logical_and(xwidth > self.min_pix, ywidth > self.min_pix),
                np.logical_and(levels >= min_level, levels <= max_level),
            )

            # Grids can either be set by edgecolors OR a colormap.
            if self.edgecolors is not None:
                edgecolors = colorConverter.to_rgba(self.edgecolors, alpha=self.alpha)
            else:  # use colormap if not explicitly overridden by edgecolors
                if self.cmap is not None:
                    color_bounds = [0, plot.data.ds.index.max_level]
                    edgecolors = (
                        apply_colormap(
                            levels[visible] * 1.0,
                            color_bounds=color_bounds,
                            cmap_name=self.cmap,
                        )[0, :, :]
                        * 1.0
                        / 255.0
                    )
                    edgecolors[:, 3] = self.alpha
                else:
                    edgecolors = (0.0, 0.0, 0.0, self.alpha)

            if visible.nonzero()[0].size == 0:
                continue

            edge_x = (left_edge_x, left_edge_x, right_edge_x, right_edge_x)
            edge_y = (left_edge_y, right_edge_y, right_edge_y, left_edge_y)
            edge_x, edge_y = self._sanitize_xy_order(plot, edge_x, edge_y)
            verts = np.array([edge_x, edge_y])
            verts = verts.transpose()[visible, :, :]
            grid_collection = matplotlib.collections.PolyCollection(
                verts,
                facecolors="none",
                edgecolors=edgecolors,
                linewidth=self.linewidth,
            )
            plot._axes.add_collection(grid_collection)

            visible_ids = np.logical_and(
                np.logical_and(xwidth > self.min_pix_ids, ywidth > self.min_pix_ids),
                np.logical_and(levels >= min_level, levels <= max_level),
            )

            if self.draw_ids:
                plot_ids = np.where(visible_ids)[0]
                x = np.empty(plot_ids.size)
                y = np.empty(plot_ids.size)
                for i, n in enumerate(plot_ids):
                    if self.id_loc == "lower left":
                        x[i] = left_edge_x[n] + (2 * (xx1 - xx0) / xpix)
                        y[i] = left_edge_y[n] + (2 * (yy1 - yy0) / ypix)
                    elif self.id_loc == "lower right":
                        x[i] = right_edge_x[n] - (
                            (10 * len(str(block_ids[i])) - 2) * (xx1 - xx0) / xpix
                        )
                        y[i] = left_edge_y[n] + (2 * (yy1 - yy0) / ypix)
                    elif self.id_loc == "upper left":
                        x[i] = left_edge_x[n] + (2 * (xx1 - xx0) / xpix)
                        y[i] = right_edge_y[n] - (12 * (yy1 - yy0) / ypix)
                    elif self.id_loc == "upper right":
                        x[i] = right_edge_x[n] - (
                            (10 * len(str(block_ids[i])) - 2) * (xx1 - xx0) / xpix
                        )
                        y[i] = right_edge_y[n] - (12 * (yy1 - yy0) / ypix)
                    else:
                        raise RuntimeError(
                            "Unrecognized id_loc value ('%s'). "
                            "Allowed values are 'lower left', lower right', "
                            "'upper left', and 'upper right'." % self.id_loc
                        )
                    xi, yi = self._sanitize_xy_order(plot, x[i], y[i])
                    plot._axes.text(xi, yi, "%d" % block_ids[n], clip_on=True)


class StreamlineCallback(PlotCallback):
    """
    Add streamlines to any plot, using the *field_x* and *field_y*
    from the associated data, skipping every *factor* datapoints like
    'quiver'. *density* is the index of the amount of the streamlines.
    *field_color* is a field to be used to colormap the streamlines.
    If *display_threshold* is supplied, any streamline segments where
    *field_color* is less than the threshold will be removed by having
    their line width set to 0.
    """

    _type_name = "streamlines"
    _supported_geometries = ("cartesian", "spectral_cube", "polar", "cylindrical")
    _incompatible_plot_types = ("OffAxisProjection", "Particle")

    def __init__(
        self,
        field_x,
        field_y,
        *,
        factor: Union[Tuple[int, int], int] = 16,
        density=1,
        field_color=None,
        display_threshold=None,
        plot_args=None,
        **kwargs,
    ):
        self.field_x = field_x
        self.field_y = field_y
        self.field_color = field_color
        self.factor = _validate_factor_tuple(factor)
        self.dens = density
        self.display_threshold = display_threshold

        if plot_args is not None:
            issue_deprecation_warning(
                "`plot_args` is deprecated. "
                "You can now pass arbitrary keyword arguments instead of a dictionary.",
                since="4.1.0",
            )
            plot_args.update(kwargs)
        else:
            plot_args = kwargs

        self.plot_args = plot_args

    def __call__(self, plot):
        bounds = self._physical_bounds(plot)
        xx0, xx1, yy0, yy1 = self._plot_bounds(plot)

        # We are feeding this size into the pixelizer, where it will properly
        # set it in reverse order
        nx = plot.raw_image_shape[1] // self.factor[0]
        ny = plot.raw_image_shape[0] // self.factor[1]
        pixX = plot.data.ds.coordinates.pixelize(
            plot.data.axis, plot.data, self.field_x, bounds, (nx, ny)
        )
        pixY = plot.data.ds.coordinates.pixelize(
            plot.data.axis, plot.data, self.field_y, bounds, (nx, ny)
        )
        if self.field_color:
            field_colors = plot.data.ds.coordinates.pixelize(
                plot.data.axis, plot.data, self.field_color, bounds, (nx, ny)
            )

            if self.display_threshold:

                mask = field_colors > self.display_threshold
                lwdefault = matplotlib.rcParams["lines.linewidth"]

                if "linewidth" in self.plot_args:
                    linewidth = self.plot_args["linewidth"]
                else:
                    linewidth = lwdefault

                try:
                    linewidth *= mask
                    self.plot_args["linewidth"] = linewidth
                except ValueError as e:
                    err_msg = (
                        "Error applying display threshold: linewidth"
                        + "must have shape ({}, {}) or be scalar"
                    )
                    err_msg = err_msg.format(nx, ny)
                    raise ValueError(err_msg) from e

        else:
            field_colors = None

        X, Y = (
            np.linspace(xx0, xx1, nx, endpoint=True),
            np.linspace(yy0, yy1, ny, endpoint=True),
        )
        X, Y, pixX, pixY = self._sanitize_xy_order(plot, X, Y, pixX, pixY)
        if plot._swap_axes:
            # need an additional transpose here for streamline tracing
            pixX = pixX.transpose()
            pixY = pixY.transpose()
            X = X.transpose()
            Y = Y.transpose()
        streamplot_args = {
            "x": X,
            "y": Y,
            "u": pixX,
            "v": pixY,
            "density": self.dens,
            "color": field_colors,
        }
        streamplot_args.update(self.plot_args)
        plot._axes.streamplot(**streamplot_args)
        self._set_plot_limits(plot, (xx0, xx1, yy0, yy1))


class LinePlotCallback(PlotCallback):
    """
    Overplot a line with endpoints at p1 and p2.  p1 and p2
    should be 2D or 3D coordinates consistent with the coordinate
    system denoted in the "coord_system" keyword.

    Parameters
    ----------
    p1, p2 : 2- or 3-element tuples, lists, or arrays
        These are the coordinates of the endpoints of the line.

    coord_system : string, optional
        This string defines the coordinate system of the coordinates p1 and p2.
        Valid coordinates are:

            "data" -- the 3D dataset coordinates

            "plot" -- the 2D coordinates defined by the actual plot limits

            "axis" -- the MPL axis coordinates: (0,0) is lower left; (1,1) is
                      upper right

            "figure" -- the MPL figure coordinates: (0,0) is lower left, (1,1)
                        is upper right

    plot_args : dictionary, optional
        This dictionary is passed to the MPL plot function for generating
        the line.  By default, it is: {'color':'white', 'linewidth':2}

    Examples
    --------

    >>> # Overplot a diagonal white line from the lower left corner to upper
    >>> # right corner
    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_line([0, 0], [1, 1], coord_system="axis")
    >>> s.save()

    >>> # Overplot a red dashed line from data coordinate (0.1, 0.2, 0.3) to
    >>> # (0.5, 0.6, 0.7)
    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_line(
    ...     [0.1, 0.2, 0.3],
    ...     [0.5, 0.6, 0.7],
    ...     coord_system="data",
    ...     color="red",
    ...     linestyles="--",
    ... )
    >>> s.save()

    """

    _type_name = "line"
    _supported_geometries: Tuple[str, ...] = (
        "cartesian",
        "spectral_cube",
        "polar",
        "cylindrical",
    )

    def __init__(
        self,
        p1,
        p2,
        *,
        coord_system="data",
        data_coords=False,
        plot_args: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        def_plot_args = {"color": "white", "linewidth": 2}

        self.p1 = p1
        self.p2 = p2
        if plot_args is not None:
            issue_deprecation_warning(
                "`plot_args` is deprecated. "
                "You can now pass arbitrary keyword arguments instead of a dictionary.",
                since="4.1.0",
            )
            plot_args.update(kwargs)

        self.plot_args = {**def_plot_args, **kwargs}

        if data_coords:
            coord_system = "data"
            issue_deprecation_warning(
                "The data_coords keyword is deprecated.  Please set "
                "the keyword coord_system='data' instead.",
                since="3.2",
                removal="4.2",
            )
        self.coord_system = coord_system
        self.transform = None

    def __call__(self, plot):
        p1 = self._sanitize_coord_system(plot, self.p1, coord_system=self.coord_system)
        p2 = self._sanitize_coord_system(plot, self.p2, coord_system=self.coord_system)
        start_pt, end_pt = [p1[0], p2[0]], [p1[1], p2[1]]
        if plot._swap_axes:
            start_pt, end_pt = [p2[0], p1[0]], [p2[1], p1[1]]
        plot._axes.plot(start_pt, end_pt, transform=self.transform, **self.plot_args)
        self._set_plot_limits(plot)


class ImageLineCallback(LinePlotCallback):
    """
    This callback is deprecated, as it is simply a wrapper around
    the LinePlotCallback (ie annotate_image()).  The only difference is
    that it uses coord_system="axis" by default. Please see LinePlotCallback
    for more information.

    """

    _type_name = "image_line"
    _supported_geometries = ("cartesian", "spectral_cube", "cylindrical")

    def __init__(self, p1, p2, data_coords=False, coord_system="axis", plot_args=None):
        super().__init__(
            p1,
            p2,
            data_coords=data_coords,
            coord_system=coord_system,
            plot_args=plot_args,
        )
        issue_deprecation_warning(
            "The ImageLineCallback (annotate_image_line()) is "
            "deprecated.  Please use the LinePlotCallback "
            "(annotate_line()) instead.",
            since="3.2",
            removal="4.2",
        )

    def __call__(self, plot):
        super().__call__(plot)


class CuttingQuiverCallback(BaseQuiverCallback):
    """
    Get a quiver plot on top of a cutting plane, using *field_x* and
    *field_y*, skipping every *factor* datapoint in the discretization.
    *scale* is the data units per arrow length unit using *scale_units*
    and *plot_args* allows you to pass in matplotlib arguments (see
    matplotlib.axes.Axes.quiver for more info). if *normalize* is True,
    the fields will be scaled by their local (in-plane) length, allowing
    morphological features to be more clearly seen for fields with
    substantial variation in field strength.
    """

    _type_name = "cquiver"
    _supported_geometries = ("cartesian", "spectral_cube")

    def _get_quiver_data(self, plot, bounds: tuple, nx: int, ny: int):
        # calls the pixelizer, returns pixX, pixY, pixC arrays
        indices = np.argsort(plot.data["index", "dx"])[::-1].astype(np.int_)

        pixX = np.zeros((ny, nx), dtype="f8")
        pixY = np.zeros((ny, nx), dtype="f8")
        pixelize_off_axis_cartesian(
            pixX,
            plot.data[("index", "x")].to("code_length"),
            plot.data[("index", "y")].to("code_length"),
            plot.data[("index", "z")].to("code_length"),
            plot.data["px"],
            plot.data["py"],
            plot.data["pdx"],
            plot.data["pdy"],
            plot.data["pdz"],
            plot.data.center,
            plot.data._inv_mat,
            indices,
            plot.data[self.field_x],
            bounds,
        )
        pixelize_off_axis_cartesian(
            pixY,
            plot.data[("index", "x")].to("code_length"),
            plot.data[("index", "y")].to("code_length"),
            plot.data[("index", "z")].to("code_length"),
            plot.data["px"],
            plot.data["py"],
            plot.data["pdx"],
            plot.data["pdy"],
            plot.data["pdz"],
            plot.data.center,
            plot.data._inv_mat,
            indices,
            plot.data[self.field_y],
            bounds,
        )

        if self.field_c is None:
            pixC = None
        else:
            pixC = np.zeros((ny, nx), dtype="f8")
            pixelize_off_axis_cartesian(
                pixC,
                plot.data[("index", "x")].to("code_length"),
                plot.data[("index", "y")].to("code_length"),
                plot.data[("index", "z")].to("code_length"),
                plot.data["px"],
                plot.data["py"],
                plot.data["pdx"],
                plot.data["pdy"],
                plot.data["pdz"],
                plot.data.center,
                plot.data._inv_mat,
                indices,
                plot.data[self.field_c],
                bounds,
            )

        return pixX, pixY, pixC


class ClumpContourCallback(PlotCallback):
    """
    Take a list of *clumps* and plot them as a set of contours.
    """

    _type_name = "clumps"
    _supported_geometries = ("cartesian", "spectral_cube", "cylindrical")
    _incompatible_plot_types = ("OffAxisSlice", "OffAxisProjection", "Particle")

    def __init__(self, clumps, *, plot_args=None, **kwargs):
        self.clumps = clumps
        if plot_args is not None:
            issue_deprecation_warning(
                "`plot_args` is deprecated. "
                "You can now pass arbitrary keyword arguments instead of a dictionary.",
                since="4.1.0",
            )
            plot_args.update(kwargs)
        else:
            plot_args = kwargs
        if "color" in plot_args:
            plot_args["colors"] = plot_args.pop("color")
        self.plot_args = plot_args

    def __call__(self, plot):
        bounds = self._physical_bounds(plot)
        extent = self._plot_bounds(plot)

        ax = plot.data.axis
        px_index = plot.data.ds.coordinates.x_axis[ax]
        py_index = plot.data.ds.coordinates.y_axis[ax]

        xf = plot.data.ds.coordinates.axis_name[px_index]
        yf = plot.data.ds.coordinates.axis_name[py_index]
        dxf = f"d{xf}"
        dyf = f"d{yf}"

        ny, nx = plot.raw_image_shape
        buff = np.zeros((nx, ny), dtype="float64")
        for i, clump in enumerate(reversed(self.clumps)):
            mylog.info("Pixelizing contour %s", i)

            if isinstance(clump, Clump):
                ftype = "index"
            elif isinstance(clump, YTClumpContainer):
                ftype = "grid"
            else:
                raise RuntimeError(
                    f"Unknown field type for object of type {type(clump)}."
                )

            xf_copy = clump[ftype, xf].copy().in_units("code_length")
            yf_copy = clump[ftype, yf].copy().in_units("code_length")

            temp = np.zeros((ny, nx), dtype="f8")
            pixelize_cartesian(
                temp,
                xf_copy,
                yf_copy,
                clump[ftype, dxf].in_units("code_length") / 2.0,
                clump[ftype, dyf].in_units("code_length") / 2.0,
                clump[ftype, dxf].d * 0.0 + i + 1,  # inits inside Pixelize
                bounds,
                0,
            )
            buff = np.maximum(temp, buff)
        if plot._swap_axes:
            buff = buff.transpose()
            extent = (extent[2], extent[3], extent[0], extent[1])
        self.rv = plot._axes.contour(
            buff, np.unique(buff), extent=extent, **self.plot_args
        )


class ArrowCallback(PlotCallback):
    """
    Overplot arrow(s) pointing at position(s) for highlighting specific
    features.  By default, arrow points from lower left to the designated
    position "pos" with arrow length "length".  Alternatively, if
    "starting_pos" is set, arrow will stretch from "starting_pos" to "pos"
    and "length" will be disregarded.

    "coord_system" keyword refers to positions set in "pos" arg and
    "starting_pos" keyword, which by default are in data coordinates.

    "length", "width", "head_length", and "head_width" keywords for the arrow
    are all in axis units, ie relative to the size of the plot axes as 1,
    even if the position of the arrow is set relative to another coordinate
    system.

    Parameters
    ----------
    pos : array-like
        These are the coordinates where the marker(s) will be overplotted
        Either as [x,y,z] or as [[x1,x2,...],[y1,y2,...],[z1,z2,...]]

    length : float, optional
        The length, in axis units, of the arrow.
        Default: 0.03

    width : float, optional
        The width, in axis units, of the tail line of the arrow.
        Default: 0.003

    head_length : float, optional
        The length, in axis units, of the head of the arrow.  If set
        to None, use 1.5*head_width
        Default: None

    head_width : float, optional
        The width, in axis units, of the head of the arrow.
        Default: 0.02

    starting_pos : 2- or 3-element tuple, list, or array, optional
        These are the coordinates from which the arrow starts towards its
        point.  Not compatible with 'length' kwarg.

    coord_system : string, optional
        This string defines the coordinate system of the coordinates of pos
        Valid coordinates are:

            "data" -- the 3D dataset coordinates

            "plot" -- the 2D coordinates defined by the actual plot limits

            "axis" -- the MPL axis coordinates: (0,0) is lower left; (1,1) is
                      upper right

            "figure" -- the MPL figure coordinates: (0,0) is lower left, (1,1)
                        is upper right

    plot_args : dictionary, optional
        This dictionary is passed to the MPL arrow function for generating
        the arrow.  By default, it is: {'color':'white'}

    Examples
    --------

    >>> # Overplot an arrow pointing to feature at data coord: (0.2, 0.3, 0.4)
    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_arrow([0.2, 0.3, 0.4])
    >>> s.save()

    >>> # Overplot a red arrow with longer length pointing to plot coordinate
    >>> # (0.1, -0.1)
    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_arrow(
    ...     [0.1, -0.1], length=0.06, coord_system="plot", color="red"
    ... )
    >>> s.save()

    """

    _type_name = "arrow"
    _supported_geometries = ("cartesian", "spectral_cube", "cylindrical")

    def __init__(
        self,
        pos,
        *,
        code_size=None,
        length=0.03,
        width=0.0001,
        head_width=0.01,
        head_length=0.01,
        starting_pos=None,
        coord_system="data",
        plot_args: Optional[Dict[str, Any]] = None,  # deprecated
        **kwargs,
    ):
        def_plot_args = {"color": "white"}
        self.pos = pos
        self.code_size = code_size
        self.length = length
        self.width = width
        self.head_width = head_width
        self.head_length = head_length
        self.starting_pos = starting_pos
        self.coord_system = coord_system
        self.transform = None

        if plot_args is not None:
            issue_deprecation_warning(
                "`plot_args` is deprecated. "
                "You can now pass arbitrary keyword arguments instead of a dictionary.",
                since="4.1.0",
            )
            plot_args = {**def_plot_args, **plot_args, **kwargs}
        else:
            plot_args = def_plot_args
        self.plot_args = plot_args

    def __call__(self, plot):
        x, y = self._sanitize_coord_system(
            plot, self.pos, coord_system=self.coord_system
        )
        xx0, xx1, yy0, yy1 = self._plot_bounds(plot)
        # normalize all of the kwarg lengths to the plot size
        plot_diag = ((yy1 - yy0) ** 2 + (xx1 - xx0) ** 2) ** (0.5)
        self.length *= plot_diag
        self.width *= plot_diag
        self.head_width *= plot_diag
        if self.head_length is not None:
            self.head_length *= plot_diag
        if self.code_size is not None:
            issue_deprecation_warning(
                "The code_size keyword is deprecated.  Please use "
                "the length keyword in 'axis' units instead. "
                "Setting code_size overrides length value.",
                since="3.2",
                removal="4.2",
            )
            if is_sequence(self.code_size):
                self.code_size = plot.data.ds.quan(self.code_size[0], self.code_size[1])
                self.code_size = np.float64(self.code_size.in_units(plot.xlim[0].units))
            self.code_size = self.code_size * self._pixel_scale(plot)[0]
            dx = dy = self.code_size
        else:
            if self.starting_pos is not None:
                start_x, start_y = self._sanitize_coord_system(
                    plot, self.starting_pos, coord_system=self.coord_system
                )
                dx = x - start_x
                dy = y - start_y
            else:
                dx = (xx1 - xx0) * 2 ** (0.5) * self.length
                dy = (yy1 - yy0) * 2 ** (0.5) * self.length
        # If the arrow is 0 length
        if dx == dy == 0:
            warnings.warn("The arrow has zero length.  Not annotating.")
            return

        x, y, dx, dy = self._sanitize_xy_order(plot, x, y, dx, dy)
        try:
            plot._axes.arrow(
                x - dx,
                y - dy,
                dx,
                dy,
                width=self.width,
                head_width=self.head_width,
                head_length=self.head_length,
                transform=self.transform,
                length_includes_head=True,
                **self.plot_args,
            )
        except ValueError:
            for i in range(len(x)):
                plot._axes.arrow(
                    x[i] - dx,
                    y[i] - dy,
                    dx,
                    dy,
                    width=self.width,
                    head_width=self.head_width,
                    head_length=self.head_length,
                    transform=self.transform,
                    length_includes_head=True,
                    **self.plot_args,
                )
        self._set_plot_limits(plot, (xx0, xx1, yy0, yy1))


class MarkerAnnotateCallback(PlotCallback):
    """
    Overplot marker(s) at a position(s) for highlighting specific features.

    Parameters
    ----------
    pos : array-like
        These are the coordinates where the marker(s) will be overplotted
        Either as [x,y,z] or as [[x1,x2,...],[y1,y2,...],[z1,z2,...]]

    marker : string, optional
        The shape of the marker to be passed to the MPL scatter function.
        By default, it is 'x', but other acceptable values are: '.', 'o', 'v',
        '^', 's', 'p' '*', etc.  See matplotlib.markers for more information.

    coord_system : string, optional
        This string defines the coordinate system of the coordinates of pos
        Valid coordinates are:

            "data" -- the 3D dataset coordinates

            "plot" -- the 2D coordinates defined by the actual plot limits

            "axis" -- the MPL axis coordinates: (0,0) is lower left; (1,1) is
                      upper right

            "figure" -- the MPL figure coordinates: (0,0) is lower left, (1,1)
                        is upper right

    plot_args : dictionary, optional
        This dictionary is passed to the MPL scatter function for generating
        the marker.  By default, it is: {'color':'white', 's':50}

    Examples
    --------

    >>> # Overplot a white X on a feature at data location (0.5, 0.5, 0.5)
    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_marker([0.4, 0.5, 0.6])
    >>> s.save()

    >>> # Overplot a big yellow circle at axis location (0.1, 0.2)
    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_marker(
    ...     [0.1, 0.2],
    ...     marker="o",
    ...     coord_system="axis",
    ...     color="yellow",
    ...     s=200,
    ... )
    >>> s.save()

    """

    _type_name = "marker"
    _supported_geometries = ("cartesian", "spectral_cube", "polar", "cylindrical")

    def __init__(
        self, pos, marker="x", *, coord_system="data", plot_args=None, **kwargs
    ):
        def_plot_args = {"color": "w", "s": 50}
        self.pos = pos
        self.marker = marker
        if plot_args is not None:
            issue_deprecation_warning(
                "`plot_args` is deprecated. "
                "You can now pass arbitrary keyword arguments instead of a dictionary.",
                since="4.1.0",
            )
            plot_args = {**def_plot_args, **plot_args, **kwargs}
        else:
            plot_args = {**def_plot_args, **kwargs}
        self.plot_args = plot_args
        self.coord_system = coord_system
        self.transform = None

    def __call__(self, plot):
        x, y = self._sanitize_coord_system(
            plot, self.pos, coord_system=self.coord_system
        )
        x, y = self._sanitize_xy_order(plot, x, y)
        plot._axes.scatter(
            x, y, marker=self.marker, transform=self.transform, **self.plot_args
        )
        self._set_plot_limits(plot)


class SphereCallback(PlotCallback):
    """
    Overplot a circle with designated center and radius with optional text.

    Parameters
    ----------
    center : 2- or 3-element tuple, list, or array
        These are the coordinates where the circle will be overplotted

    radius : YTArray, float, or (1, ('kpc')) style tuple
        The radius of the circle in code coordinates

    circle_args : dict, optional
        This dictionary is passed to the MPL circle object. By default,
        {'color':'white'}

    coord_system : string, optional
        This string defines the coordinate system of the coordinates of pos
        Valid coordinates are:

            "data" -- the 3D dataset coordinates

            "plot" -- the 2D coordinates defined by the actual plot limits

            "axis" -- the MPL axis coordinates: (0,0) is lower left; (1,1) is
                      upper right

            "figure" -- the MPL figure coordinates: (0,0) is lower left, (1,1)
                        is upper right

    text : string, optional
        Optional text to include next to the circle.

    text_args : dictionary, optional
        This dictionary is passed to the MPL text function. By default,
        it is: {'color':'white'}

    Examples
    --------

    >>> # Overplot a white circle of radius 100 kpc over the central galaxy
    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_sphere([0.5, 0.5, 0.5], radius=(100, "kpc"))
    >>> s.save()

    """

    _type_name = "sphere"
    _supported_geometries = ("cartesian", "spectral_cube", "polar", "cylindrical")

    def __init__(
        self,
        center,
        radius,
        *,
        coord_system="data",
        text=None,
        circle_args=None,
        text_args=None,
    ):
        def_text_args = {"color": "white"}
        def_circle_args = {"color": "white"}
        self.center = center
        self.radius = radius
        if circle_args is None:
            circle_args = def_circle_args
        if "fill" not in circle_args:
            circle_args["fill"] = False
        self.circle_args = circle_args
        self.text = text
        if text_args is None:
            text_args = def_text_args
        self.text_args = text_args
        self.coord_system = coord_system
        self.transform = None

    def __call__(self, plot):
        from matplotlib.patches import Circle

        if is_sequence(self.radius):
            self.radius = plot.data.ds.quan(self.radius[0], self.radius[1])
            self.radius = np.float64(self.radius.in_units(plot.xlim[0].units))
        if isinstance(self.radius, YTQuantity):
            if isinstance(self.center, YTArray):
                units = self.center.units
            else:
                units = "code_length"
            self.radius = self.radius.to(units)

        # This assures the radius has the appropriate size in
        # the different coordinate systems, since one cannot simply
        # apply a different transform for a length in the same way
        # you can for a coordinate.
        if self.coord_system == "data" or self.coord_system == "plot":
            self.radius = self.radius * self._pixel_scale(plot)[0]
        else:
            self.radius /= (plot.xlim[1] - plot.xlim[0]).v

        x, y = self._sanitize_coord_system(
            plot, self.center, coord_system=self.coord_system
        )

        x, y = self._sanitize_xy_order(plot, x, y)
        cir = Circle((x, y), self.radius, transform=self.transform, **self.circle_args)

        plot._axes.add_patch(cir)
        if self.text is not None:
            label = plot._axes.text(
                x, y, self.text, transform=self.transform, **self.text_args
            )
            self._set_font_properties(plot, [label], **self.text_args)

        self._set_plot_limits(plot)


class TextLabelCallback(PlotCallback):
    """
    Overplot text on the plot at a specified position. If you desire an inset
    box around your text, set one with the inset_box_args dictionary
    keyword.

    Parameters
    ----------
    pos : 2- or 3-element tuple, list, or array
        These are the coordinates where the text will be overplotted

    text : string
        The text you wish to include

    coord_system : string, optional
        This string defines the coordinate system of the coordinates of pos
        Valid coordinates are:

            "data" -- the 3D dataset coordinates

            "plot" -- the 2D coordinates defined by the actual plot limits

            "axis" -- the MPL axis coordinates: (0,0) is lower left; (1,1) is
                      upper right

            "figure" -- the MPL figure coordinates: (0,0) is lower left, (1,1)
                        is upper right

    text_args : dictionary, optional
        This dictionary is passed to the MPL text function for generating
        the text.  By default, it is: {'color':'white'} and uses the defaults
        for the other fonts in the image.

    inset_box_args : dictionary, optional
        A dictionary of any arbitrary parameters to be passed to the Matplotlib
        FancyBboxPatch object as the inset box around the text.  Default: {}

    Examples
    --------

    >>> # Overplot white text at data location [0.55, 0.7, 0.4]
    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_text([0.55, 0.7, 0.4], "Here is a galaxy")
    >>> s.save()

    >>> # Overplot yellow text at axis location [0.2, 0.8] with
    >>> # a shaded inset box
    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_text(
    ...     [0.2, 0.8],
    ...     "Here is a galaxy",
    ...     coord_system="axis",
    ...     text_args={"color": "yellow"},
    ...     inset_box_args={
    ...         "boxstyle": "square,pad=0.3",
    ...         "facecolor": "black",
    ...         "linewidth": 3,
    ...         "edgecolor": "white",
    ...         "alpha": 0.5,
    ...     },
    ... )
    >>> s.save()
    """

    _type_name = "text"
    _supported_geometries: Tuple[str, ...] = (
        "cartesian",
        "spectral_cube",
        "polar",
        "cylindrical",
    )

    def __init__(
        self,
        pos,
        text,
        *,
        coord_system="data",
        text_args=None,
        inset_box_args=None,
        data_coords=False,  # deprecated
    ):
        def_text_args = {"color": "white"}
        self.pos = pos
        self.text = text
        if data_coords:
            coord_system = "data"
            issue_deprecation_warning(
                "The data_coords keyword is deprecated.  Please set "
                "the keyword coord_system='data' instead.",
                since="3.2",
                removal="4.2",
            )
        if text_args is None:
            text_args = def_text_args
        self.text_args = text_args
        self.inset_box_args = inset_box_args
        self.coord_system = coord_system
        self.transform = None

    def __call__(self, plot):
        kwargs = self.text_args.copy()
        x, y = self._sanitize_coord_system(
            plot, self.pos, coord_system=self.coord_system
        )

        # Set the font properties of text from this callback to be
        # consistent with other text labels in this figure
        if self.inset_box_args is not None:
            kwargs["bbox"] = self.inset_box_args
        x, y = self._sanitize_xy_order(plot, x, y)
        label = plot._axes.text(x, y, self.text, transform=self.transform, **kwargs)
        self._set_font_properties(plot, [label], **kwargs)
        self._set_plot_limits(plot)


class PointAnnotateCallback(TextLabelCallback):
    """
    This callback is deprecated, as it is simply a wrapper around
    the TextLabelCallback (ie annotate_text()).  Please see TextLabelCallback
    for more information.

    """

    _type_name = "point"
    _supported_geometries = ("cartesian", "spectral_cube", "cylindrical")

    def __init__(
        self,
        pos,
        text,
        *,
        coord_system="data",
        text_args=None,
        inset_box_args=None,
        data_coords=False,  # deprecated
    ):
        super().__init__(
            pos, text, data_coords, coord_system, text_args, inset_box_args
        )
        issue_deprecation_warning(
            "The PointAnnotateCallback (annotate_point()) is "
            "deprecated.  Please use the TextLabelCallback "
            "(annotate_point()) instead.",
            since="3.2",
            removal="4.2",
        )

    def __call__(self, plot):
        super().__call__(plot)


class HaloCatalogCallback(PlotCallback):
    """
    Plots circles at the locations of all the halos
    in a halo catalog with radii corresponding to the
    virial radius of each halo.

    Note, this functionality requires the yt_astro_analysis
    package. See https://yt-astro-analysis.readthedocs.io/
    for more information.

    Parameters
    ----------
    halo_catalog : Dataset, DataContainer,
                   or ~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog
        The object containing halos to be overplotted. This can
        be a HaloCatalog object, a loaded halo catalog dataset,
        or a data container from a halo catalog dataset.
    circle_args : list
        Contains the arguments controlling the
        appearance of the circles, supplied to the
        Matplotlib patch Circle.
    width : tuple
        The width over which to select halos to plot,
        useful when overplotting to a slice plot. Accepts
        a tuple in the form (1.0, 'Mpc').
    annotate_field : str
        A field contained in the
        halo catalog to add text to the plot near the halo.
        Example: annotate_field = 'particle_mass' will
        write the halo mass next to each halo.
    radius_field : str
        A field contained in the halo
        catalog to set the radius of the circle which will
        surround each halo. Default: 'virial_radius'.
    center_field_prefix : str
        Accepts a field prefix which will
        be used to find the fields containing the coordinates
        of the center of each halo. Ex: 'particle_position'
        will result in the fields 'particle_position_x' for x
        'particle_position_y' for y, and 'particle_position_z'
        for z. Default: 'particle_position'.
    text_args : dict
        Contains the arguments controlling the text
        appearance of the annotated field.
    factor : float
        A number the virial radius is multiplied by for
        plotting the circles. Ex: factor = 2.0 will plot
        circles with twice the radius of each halo virial radius.

    Examples
    --------

    >>> import yt
    >>> dds = yt.load("Enzo_64/DD0043/data0043")
    >>> hds = yt.load("rockstar_halos/halos_0.0.bin")
    >>> p = yt.ProjectionPlot(
    ...     dds, "x", ("gas", "density"), weight_field=("gas", "density")
    ... )
    >>> p.annotate_halos(hds)
    >>> p.save()

    >>> # plot a subset of all halos
    >>> import yt
    >>> dds = yt.load("Enzo_64/DD0043/data0043")
    >>> hds = yt.load("rockstar_halos/halos_0.0.bin")
    >>> # make a region half the width of the box
    >>> dregion = dds.box(
    ...     dds.domain_center - 0.25 * dds.domain_width,
    ...     dds.domain_center + 0.25 * dds.domain_width,
    ... )
    >>> hregion = hds.box(
    ...     hds.domain_center - 0.25 * hds.domain_width,
    ...     hds.domain_center + 0.25 * hds.domain_width,
    ... )
    >>> p = yt.ProjectionPlot(
    ...     dds,
    ...     "x",
    ...     ("gas", "density"),
    ...     weight_field=("gas", "density"),
    ...     data_source=dregion,
    ...     width=0.5,
    ... )
    >>> p.annotate_halos(hregion)
    >>> p.save()

    >>> # plot halos from a HaloCatalog
    >>> import yt
    >>> from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog
    >>> dds = yt.load("Enzo_64/DD0043/data0043")
    >>> hds = yt.load("rockstar_halos/halos_0.0.bin")
    >>> hc = HaloCatalog(data_ds=dds, halos_ds=hds)
    >>> p = yt.ProjectionPlot(
    ...     dds, "x", ("gas", "density"), weight_field=("gas", "density")
    ... )
    >>> p.annotate_halos(hc)
    >>> p.save()

    """

    _type_name = "halos"
    region = None
    _descriptor = None
    _supported_geometries = ("cartesian", "spectral_cube")

    def __init__(
        self,
        halo_catalog,
        circle_args=None,
        circle_kwargs=None,
        width=None,
        annotate_field=None,
        radius_field="virial_radius",
        center_field_prefix="particle_position",
        text_args=None,
        font_kwargs=None,
        factor=1.0,
    ):
        issue_deprecation_warning(
            "The annotate_halos method has been fully migrated to the "
            "yt_astro_analysis extension. "
            "Please update the extension to version 1.1 or newer. "
            "This duplicated functionality will be removed from the main yt package.",
            since="4.1",
            removal="4.2",
        )

        try:
            from yt_astro_analysis.halo_analysis.api import HaloCatalog
        except ImportError:
            HaloCatalog = NotAModule("yt_astro_analysis")

        def_circle_args = {"edgecolor": "white", "facecolor": "None"}
        def_text_args = {"color": "white"}

        if isinstance(halo_catalog, YTDataContainer):
            self.halo_data = halo_catalog
        elif isinstance(halo_catalog, Dataset):
            self.halo_data = halo_catalog.all_data()
        elif isinstance(halo_catalog, HaloCatalog):
            if halo_catalog.data_source.ds == halo_catalog.halos_ds:
                self.halo_data = halo_catalog.data_source
            else:
                self.halo_data = halo_catalog.halos_ds.all_data()
        else:
            raise RuntimeError(
                "halo_catalog argument must be a HaloCatalog object, "
                + "a dataset, or a data container."
            )

        self.width = width
        self.radius_field = radius_field
        self.center_field_prefix = center_field_prefix
        self.annotate_field = annotate_field
        if circle_kwargs is not None:
            circle_args = circle_kwargs
            issue_deprecation_warning(
                "The circle_kwargs keyword is deprecated.  Please "
                "use the circle_args keyword instead.",
                since="3.2",
                removal="4.2",
            )
        if font_kwargs is not None:
            text_args = font_kwargs
            issue_deprecation_warning(
                "The font_kwargs keyword is deprecated.  Please use "
                "the text_args keyword instead.",
                since="3.2",
                removal="4.2",
            )
        if circle_args is None:
            circle_args = def_circle_args
        self.circle_args = circle_args
        if text_args is None:
            text_args = def_text_args
        self.text_args = text_args
        self.factor = factor

    def __call__(self, plot):
        from matplotlib.patches import Circle

        data = plot.data

        halo_data = self.halo_data
        axis_names = plot.data.ds.coordinates.axis_name
        xax = plot.data.ds.coordinates.x_axis[data.axis]
        yax = plot.data.ds.coordinates.y_axis[data.axis]
        field_x = f"{self.center_field_prefix}_{axis_names[xax]}"
        field_y = f"{self.center_field_prefix}_{axis_names[yax]}"
        field_z = f"{self.center_field_prefix}_{axis_names[data.axis]}"

        # Set up scales for pixel size and original data
        pixel_scale = self._pixel_scale(plot)[0]
        units = plot.xlim[0].units

        # Convert halo positions to code units of the plotted data
        # and then to units of the plotted window
        px = halo_data[("all", field_x)][:].in_units(units)
        py = halo_data[("all", field_y)][:].in_units(units)

        xplotcenter = (plot.xlim[0] + plot.xlim[1]) / 2
        yplotcenter = (plot.ylim[0] + plot.ylim[1]) / 2

        xdomaincenter = plot.ds.domain_center[xax]
        ydomaincenter = plot.ds.domain_center[yax]

        xoffset = xplotcenter - xdomaincenter
        yoffset = yplotcenter - ydomaincenter

        xdw = plot.ds.domain_width[xax].to(units)
        ydw = plot.ds.domain_width[yax].to(units)

        modpx = np.mod(px - xoffset, xdw) + xoffset
        modpy = np.mod(py - yoffset, ydw) + yoffset

        px[modpx != px] = modpx[modpx != px]
        py[modpy != py] = modpy[modpy != py]

        px, py = self._convert_to_plot(plot, [px, py])

        # Convert halo radii to a radius in pixels
        radius = halo_data[("all", self.radius_field)][:].in_units(units)
        radius = np.array(radius * pixel_scale * self.factor)

        if self.width:
            pz = halo_data[("all", field_z)][:].in_units("code_length")
            c = data.center[data.axis]

            # I should catch an error here if width isn't in this form
            # but I dont really want to reimplement get_sanitized_width...
            width = data.ds.arr(self.width[0], self.width[1]).in_units("code_length")

            indices = np.where((pz > c - 0.5 * width) & (pz < c + 0.5 * width))

            px = px[indices]
            py = py[indices]
            radius = radius[indices]

        px, py = self._sanitize_xy_order(plot, px, py)
        for x, y, r in zip(px, py, radius):
            plot._axes.add_artist(Circle(xy=(x, y), radius=r, **self.circle_args))
        self._set_plot_limits(plot)

        if self.annotate_field:
            annotate_dat = halo_data[("all", self.annotate_field)]
            texts = [f"{float(dat):g}" for dat in annotate_dat]
            labels = []
            for pos_x, pos_y, t in zip(px, py, texts):
                labels.append(plot._axes.text(pos_x, pos_y, t, **self.text_args))

            # Set the font properties of text from this callback to be
            # consistent with other text labels in this figure
            self._set_font_properties(plot, labels, **self.text_args)


class ParticleCallback(PlotCallback):
    """
    Adds particle positions, based on a thick slab along *axis* with a
    *width* along the line of sight.  *p_size* controls the number of
    pixels per particle, and *col* governs the color.  *ptype* will
    restrict plotted particles to only those that are of a given type.
    *alpha* determines the opacity of the marker symbol used in the scatter.
    An alternate data source can be specified with *data_source*, but by
    default the plot's data source will be queried.
    """

    _type_name = "particles"
    region = None
    _descriptor = None
    _supported_geometries = ("cartesian", "spectral_cube", "cylindrical")
    _incompatible_plot_types = ("OffAxisSlice", "OffAxisProjection")

    def __init__(
        self,
        width,
        p_size=1.0,
        col="k",
        marker="o",
        stride=1,
        ptype="all",
        minimum_mass=None,
        alpha=1.0,
        data_source=None,
    ):
        self.width = width
        self.p_size = p_size
        self.color = col
        self.marker = marker
        self.stride = stride
        self.ptype = ptype
        self.minimum_mass = minimum_mass
        self.alpha = alpha
        self.data_source = data_source
        if self.minimum_mass is not None:
            issue_deprecation_warning(
                "The minimum_mass keyword is deprecated.  Please use "
                "an appropriate particle filter and the ptype keyword instead.",
                since="3.5",
                removal="4.2",
            )

    def __call__(self, plot):
        data = plot.data
        if is_sequence(self.width):
            validate_width_tuple(self.width)
            self.width = plot.data.ds.quan(self.width[0], self.width[1])
        elif isinstance(self.width, YTQuantity):
            self.width = plot.data.ds.quan(self.width.value, self.width.units)
        else:
            self.width = plot.data.ds.quan(self.width, "code_length")
        # we construct a rectangular prism
        x0, x1, y0, y1 = self._physical_bounds(plot)
        if isinstance(self.data_source, YTCutRegion):
            mylog.warning(
                "Parameter 'width' is ignored in annotate_particles if the "
                "data_source is a cut_region. "
                "See https://github.com/yt-project/yt/issues/1933 for further details."
            )
            self.region = self.data_source
        else:
            self.region = self._get_region((x0, x1), (y0, y1), plot.data.axis, data)
        ax = data.axis
        xax = plot.data.ds.coordinates.x_axis[ax]
        yax = plot.data.ds.coordinates.y_axis[ax]
        axis_names = plot.data.ds.coordinates.axis_name
        field_x = f"particle_position_{axis_names[xax]}"
        field_y = f"particle_position_{axis_names[yax]}"
        pt = self.ptype
        self.periodic_x = plot.data.ds.periodicity[xax]
        self.periodic_y = plot.data.ds.periodicity[yax]
        self.LE = plot.data.ds.domain_left_edge[xax], plot.data.ds.domain_left_edge[yax]
        self.RE = (
            plot.data.ds.domain_right_edge[xax],
            plot.data.ds.domain_right_edge[yax],
        )
        period_x = plot.data.ds.domain_width[xax]
        period_y = plot.data.ds.domain_width[yax]
        particle_x, particle_y = self._enforce_periodic(
            self.region[pt, field_x],
            self.region[pt, field_y],
            x0,
            x1,
            period_x,
            y0,
            y1,
            period_y,
        )
        gg = (
            (particle_x >= x0)
            & (particle_x <= x1)
            & (particle_y >= y0)
            & (particle_y <= y1)
        )
        if self.minimum_mass is not None:
            gg &= self.region[pt, "particle_mass"] >= self.minimum_mass
            if gg.sum() == 0:
                return
        px, py = [particle_x[gg][:: self.stride], particle_y[gg][:: self.stride]]
        px, py = self._convert_to_plot(plot, [px, py])
        px, py = self._sanitize_xy_order(plot, px, py)
        plot._axes.scatter(
            px,
            py,
            edgecolors="None",
            marker=self.marker,
            s=self.p_size,
            c=self.color,
            alpha=self.alpha,
        )
        self._set_plot_limits(plot)

    def _enforce_periodic(
        self, particle_x, particle_y, x0, x1, period_x, y0, y1, period_y
    ):
        #  duplicate particles if periodic in that direction AND if the plot
        #  extends outside the domain boundaries.
        if self.periodic_x and x0 > self.RE[0]:
            particle_x = uhstack((particle_x, particle_x + period_x))
            particle_y = uhstack((particle_y, particle_y))
        if self.periodic_x and x1 < self.LE[0]:
            particle_x = uhstack((particle_x, particle_x - period_x))
            particle_y = uhstack((particle_y, particle_y))
        if self.periodic_y and y0 > self.RE[1]:
            particle_y = uhstack((particle_y, particle_y + period_y))
            particle_x = uhstack((particle_x, particle_x))
        if self.periodic_y and y1 < self.LE[1]:
            particle_y = uhstack((particle_y, particle_y - period_y))
            particle_x = uhstack((particle_x, particle_x))
        return particle_x, particle_y

    def _get_region(self, xlim, ylim, axis, data):
        LE, RE = [None] * 3, [None] * 3
        ds = data.ds
        xax = ds.coordinates.x_axis[axis]
        yax = ds.coordinates.y_axis[axis]
        zax = axis
        LE[xax], RE[xax] = xlim
        LE[yax], RE[yax] = ylim
        LE[zax] = data.center[zax] - self.width * 0.5
        LE[zax].convert_to_units("code_length")
        RE[zax] = LE[zax] + self.width
        if (
            self.region is not None
            and np.all(self.region.left_edge <= LE)
            and np.all(self.region.right_edge >= RE)
        ):
            return self.region
        self.region = data.ds.region(data.center, LE, RE, data_source=self.data_source)
        return self.region


class TitleCallback(PlotCallback):
    """
    Accepts a *title* and adds it to the plot
    """

    _type_name = "title"

    def __init__(self, title):
        self.title = title

    def __call__(self, plot):
        plot._axes.set_title(self.title)
        # Set the font properties of text from this callback to be
        # consistent with other text labels in this figure
        label = plot._axes.title
        self._set_font_properties(plot, [label])


class MeshLinesCallback(PlotCallback):
    """
    Adds mesh lines to the plot. Only works for unstructured or
    semi-structured mesh data. For structured grid data, see
    GridBoundaryCallback or CellEdgesCallback.

    Parameters
    ----------

    plot_args:   dict, optional
        A dictionary of arguments that will be passed to matplotlib.

    Example
    -------

    >>> import yt
    >>> ds = yt.load("MOOSE_sample_data/out.e-s010")
    >>> sl = yt.SlicePlot(ds, "z", ("connect2", "convected"))
    >>> sl.annotate_mesh_lines(color="black")

    """

    _type_name = "mesh_lines"
    _supported_geometries = ("cartesian", "spectral_cube")
    _incompatible_plot_types = ("OffAxisSlice", "OffAxisProjection")

    def __init__(self, *, plot_args=None, **kwargs):
        if plot_args is not None:
            issue_deprecation_warning(
                "`plot_args` is deprecated. "
                "You can now pass arbitrary keyword arguments instead of a dictionary.",
                since="4.1.0",
            )
            plot_args.update(kwargs)
        else:
            plot_args = kwargs
        self.plot_args = plot_args

    def promote_2d_to_3d(self, coords, indices, plot):
        new_coords = np.zeros((2 * coords.shape[0], 3))
        new_connects = np.zeros(
            (indices.shape[0], 2 * indices.shape[1]), dtype=np.int64
        )

        new_coords[0 : coords.shape[0], 0:2] = coords
        new_coords[0 : coords.shape[0], 2] = plot.ds.domain_left_edge[2]
        new_coords[coords.shape[0] :, 0:2] = coords
        new_coords[coords.shape[0] :, 2] = plot.ds.domain_right_edge[2]

        new_connects[:, 0 : indices.shape[1]] = indices
        new_connects[:, indices.shape[1] :] = indices + coords.shape[0]

        return new_coords, new_connects

    def __call__(self, plot):

        index = plot.ds.index
        if not issubclass(type(index), UnstructuredIndex):
            raise RuntimeError(
                "Mesh line annotations only work for "
                "unstructured or semi-structured mesh data."
            )
        for i, m in enumerate(index.meshes):
            try:
                ftype, fname = plot.field
                if ftype.startswith("connect") and int(ftype[-1]) - 1 != i:
                    continue
            except ValueError:
                pass
            coords = m.connectivity_coords
            indices = m.connectivity_indices - m._index_offset

            num_verts = indices.shape[1]
            num_dims = coords.shape[1]

            if num_dims == 2 and num_verts == 3:
                coords, indices = self.promote_2d_to_3d(coords, indices, plot)
            elif num_dims == 2 and num_verts == 4:
                coords, indices = self.promote_2d_to_3d(coords, indices, plot)

            tri_indices = triangulate_indices(indices.astype(np.int_))
            points = coords[tri_indices]

            tfc = TriangleFacetsCallback(points, **self.plot_args)
            tfc(plot)


class TriangleFacetsCallback(PlotCallback):
    """
    Intended for representing a slice of a triangular faceted
    geometry in a slice plot.

    Uses a set of *triangle_vertices* to find all triangles the plane of a
    SlicePlot intersects with. The lines between the intersection points
    of the triangles are then added to the plot to create an outline
    of the geometry represented by the triangles.
    """

    _type_name = "triangle_facets"
    _supported_geometries = ("cartesian", "spectral_cube")

    def __init__(self, triangle_vertices, *, plot_args=None, **kwargs):
        if plot_args is not None:
            issue_deprecation_warning(
                "`plot_args` is deprecated. "
                "You can now pass arbitrary keyword arguments instead of a dictionary.",
                since="4.1.0",
            )
            plot_args.update(kwargs)
        else:
            plot_args = kwargs
        self.plot_args = kwargs
        self.vertices = triangle_vertices

    def __call__(self, plot):
        ax = plot.data.axis
        xax = plot.data.ds.coordinates.x_axis[ax]
        yax = plot.data.ds.coordinates.y_axis[ax]

        if not hasattr(self.vertices, "in_units"):
            vertices = plot.data.pf.arr(self.vertices, "code_length")
        else:
            vertices = self.vertices
        l_cy = triangle_plane_intersect(plot.data.axis, plot.data.coord, vertices)[
            :, :, (xax, yax)
        ]
        # l_cy is shape (nlines, 2, 2)
        # reformat for conversion to plot coordinates
        l_cy = np.rollaxis(l_cy, 0, 3)  # shape is now (2, 2, nlines)
        # convert all line starting points
        l_cy[0] = self._convert_to_plot(plot, l_cy[0])
        # convert all line ending points
        l_cy[1] = self._convert_to_plot(plot, l_cy[1])
        if plot._swap_axes:
            # more convenient to swap the x, y values here before final roll
            x0, y0 = l_cy[0]  # x, y values of start points
            x1, y1 = l_cy[1]  # x, y values of end points
            l_cy[0] = np.row_stack([y0, x0])  # swap x, y for start points
            l_cy[1] = np.row_stack([y1, x1])  # swap x, y for end points
        # convert back to shape (nlines, 2, 2)
        l_cy = np.rollaxis(l_cy, 2, 0)
        # create line collection and add it to the plot
        lc = matplotlib.collections.LineCollection(l_cy, **self.plot_args)
        plot._axes.add_collection(lc)


class TimestampCallback(PlotCallback):
    r"""
    Annotates the timestamp and/or redshift of the data output at a specified
    location in the image (either in a present corner, or by specifying (x,y)
    image coordinates with the x_pos, y_pos arguments.  If no time_units are
    specified, it will automatically choose appropriate units.  It allows for
    custom formatting of the time and redshift information, as well as the
    specification of an inset box around the text.

    Parameters
    ----------

    x_pos, y_pos : floats, optional
        The image location of the timestamp in the coord system defined by the
        coord_system kwarg.  Setting x_pos and y_pos overrides the corner
        parameter.

    corner : string, optional
        Corner sets up one of 4 predeterimined locations for the timestamp
        to be displayed in the image: 'upper_left', 'upper_right', 'lower_left',
        'lower_right' (also allows None). This value will be overridden by the
        optional x_pos and y_pos keywords.

    time : boolean, optional
        Whether or not to show the ds.current_time of the data output.  Can
        be used solo or in conjunction with redshift parameter.

    redshift : boolean, optional
        Whether or not to show the ds.current_time of the data output.  Can
        be used solo or in conjunction with the time parameter.

    time_format : string, optional
        This specifies the format of the time output assuming "time" is the
        number of time and "unit" is units of the time (e.g. 's', 'Myr', etc.)
        The time can be specified to arbitrary precision according to printf
        formatting codes (defaults to .1f -- a float with 1 digits after
        decimal).  Example: "Age = {time:.2f} {units}".

    time_unit : string, optional
        time_unit must be a valid yt time unit (e.g. 's', 'min', 'hr', 'yr',
        'Myr', etc.)

    redshift_format : string, optional
        This specifies the format of the redshift output.  The redshift can
        be specified to arbitrary precision according to printf formatting
        codes (defaults to 0.2f -- a float with 2 digits after decimal).
        Example: "REDSHIFT = {redshift:03.3g}",

    draw_inset_box : boolean, optional
        Whether or not an inset box should be included around the text
        If so, it uses the inset_box_args to set the matplotlib FancyBboxPatch
        object.

    coord_system : string, optional
        This string defines the coordinate system of the coordinates of pos
        Valid coordinates are:

        - "data": 3D dataset coordinates
        - "plot": 2D coordinates defined by the actual plot limits
        - "axis": MPL axis coordinates: (0,0) is lower left; (1,1) is upper right
        - "figure": MPL figure coordinates: (0,0) is lower left, (1,1) is upper right

    time_offset : float, (value, unit) tuple, or YTQuantity, optional
        Apply an offset to the time shown in the annotation from the
        value of the current time. If a scalar value with no units is
        passed in, the value of the *time_unit* kwarg is used for the
        units. Default: None, meaning no offset.

    text_args : dictionary, optional
        A dictionary of any arbitrary parameters to be passed to the Matplotlib
        text object.  Defaults: ``{'color':'white',
        'horizontalalignment':'center', 'verticalalignment':'top'}``.

    inset_box_args : dictionary, optional
        A dictionary of any arbitrary parameters to be passed to the Matplotlib
        FancyBboxPatch object as the inset box around the text.
        Defaults: ``{'boxstyle':'square', 'pad':0.3, 'facecolor':'black',
        'linewidth':3, 'edgecolor':'white', 'alpha':0.5}``

    Example
    -------

    >>> import yt
    >>> ds = yt.load("Enzo_64/DD0020/data0020")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_timestamp()
    """

    _type_name = "timestamp"
    _supported_geometries = ("cartesian", "spectral_cube", "cylindrical")

    def __init__(
        self,
        x_pos=None,
        y_pos=None,
        corner="lower_left",
        *,
        time=True,
        redshift=False,
        time_format="t = {time:.1f} {units}",
        time_unit=None,
        redshift_format="z = {redshift:.2f}",
        draw_inset_box=False,
        coord_system="axis",
        time_offset=None,
        text_args=None,
        inset_box_args=None,
    ):

        def_text_args = {
            "color": "white",
            "horizontalalignment": "center",
            "verticalalignment": "top",
        }
        def_inset_box_args = {
            "boxstyle": "square,pad=0.3",
            "facecolor": "black",
            "linewidth": 3,
            "edgecolor": "white",
            "alpha": 0.5,
        }

        # Set position based on corner argument.
        self.pos = (x_pos, y_pos)
        self.corner = corner
        self.time = time
        self.redshift = redshift
        self.time_format = time_format
        self.redshift_format = redshift_format
        self.time_unit = time_unit
        self.coord_system = coord_system
        self.time_offset = time_offset
        if text_args is None:
            text_args = def_text_args
        self.text_args = text_args
        if inset_box_args is None:
            inset_box_args = def_inset_box_args
        self.inset_box_args = inset_box_args

        # if inset box is not desired, set inset_box_args to {}
        if not draw_inset_box:
            self.inset_box_args = None

    def __call__(self, plot):
        # Setting pos overrides corner argument
        if self.pos[0] is None or self.pos[1] is None:
            if self.corner == "upper_left":
                self.pos = (0.03, 0.96)
                self.text_args["horizontalalignment"] = "left"
                self.text_args["verticalalignment"] = "top"
            elif self.corner == "upper_right":
                self.pos = (0.97, 0.96)
                self.text_args["horizontalalignment"] = "right"
                self.text_args["verticalalignment"] = "top"
            elif self.corner == "lower_left":
                self.pos = (0.03, 0.03)
                self.text_args["horizontalalignment"] = "left"
                self.text_args["verticalalignment"] = "bottom"
            elif self.corner == "lower_right":
                self.pos = (0.97, 0.03)
                self.text_args["horizontalalignment"] = "right"
                self.text_args["verticalalignment"] = "bottom"
            elif self.corner is None:
                self.pos = (0.5, 0.5)
                self.text_args["horizontalalignment"] = "center"
                self.text_args["verticalalignment"] = "center"
            else:
                raise ValueError(
                    "Argument 'corner' must be set to "
                    "'upper_left', 'upper_right', 'lower_left', "
                    "'lower_right', or None"
                )

        self.text = ""

        # If we're annotating the time, put it in the correct format
        if self.time:
            # If no time_units are set, then identify a best fit time unit
            if self.time_unit is None:
                if plot.ds._uses_code_time_unit:
                    # if the unit system is in code units
                    # we should not convert to seconds for the plot.
                    self.time_unit = plot.ds.unit_system.base_units[dimensions.time]
                else:
                    # in the case of non- code units then we
                    self.time_unit = plot.ds.get_smallest_appropriate_unit(
                        plot.ds.current_time, quantity="time"
                    )
            t = plot.ds.current_time.in_units(self.time_unit)
            if self.time_offset is not None:
                if isinstance(self.time_offset, tuple):
                    toffset = plot.ds.quan(self.time_offset[0], self.time_offset[1])
                elif isinstance(self.time_offset, Number):
                    toffset = plot.ds.quan(self.time_offset, self.time_unit)
                elif not isinstance(self.time_offset, YTQuantity):
                    raise RuntimeError(
                        "'time_offset' must be a float, tuple, or YTQuantity!"
                    )
                t -= toffset.in_units(self.time_unit)
            try:
                # here the time unit will be in brackets on the annotation.
                un = self.time_unit.latex_representation()
                time_unit = r"$\ \ (" + un + r")$"
            except AttributeError:
                time_unit = str(self.time_unit).replace("_", " ")
            self.text += self.time_format.format(time=float(t), units=time_unit)

        # If time and redshift both shown, do one on top of the other
        if self.time and self.redshift:
            self.text += "\n"

        # If we're annotating the redshift, put it in the correct format
        if self.redshift:
            try:
                z = plot.data.ds.current_redshift
            except AttributeError:
                raise AttributeError(
                    "Dataset does not have current_redshift. Set redshift=False."
                )
            # Replace instances of -0.0* with 0.0* to avoid
            # negative null redshifts (e.g., "-0.00").
            self.text += self.redshift_format.format(redshift=float(z))
            self.text = re.sub("-(0.0*)$", r"\g<1>", self.text)

        # This is just a fancy wrapper around the TextLabelCallback
        tcb = TextLabelCallback(
            self.pos,
            self.text,
            coord_system=self.coord_system,
            text_args=self.text_args,
            inset_box_args=self.inset_box_args,
        )
        return tcb(plot)


class ScaleCallback(PlotCallback):
    r"""
    Annotates the scale of the plot at a specified location in the image
    (either in a preset corner, or by specifying (x,y) image coordinates with
    the pos argument.  Coeff and units (e.g. 1 Mpc or 100 kpc) refer to the
    distance scale you desire to show on the plot.  If no coeff and units are
    specified, an appropriate pair will be determined such that your scale bar
    is never smaller than min_frac or greater than max_frac of your plottable
    axis length.  Additional customization of the scale bar is possible by
    adjusting the text_args and size_bar_args dictionaries.  The text_args
    dictionary accepts matplotlib's font_properties arguments to override
    the default font_properties for the current plot.  The size_bar_args
    dictionary accepts keyword arguments for the AnchoredSizeBar class in
    matplotlib's axes_grid toolkit.

    Parameters
    ----------

    corner : string, optional
        Corner sets up one of 4 predeterimined locations for the scale bar
        to be displayed in the image: 'upper_left', 'upper_right', 'lower_left',
        'lower_right' (also allows None). This value will be overridden by the
        optional 'pos' keyword.

    coeff : float, optional
        The coefficient of the unit defining the distance scale (e.g. 10 kpc or
        100 Mpc) for overplotting.  If set to None along with unit keyword,
        coeff will be automatically determined to be a power of 10
        relative to the best-fit unit.

    unit : string, optional
        unit must be a valid yt distance unit (e.g. 'm', 'km', 'AU', 'pc',
        'kpc', etc.) or set to None.  If set to None, will be automatically
        determined to be the best-fit to the data.

    pos : 2- or 3-element tuples, lists, or arrays, optional
        The image location of the scale bar in the plot coordinate system.
        Setting pos overrides the corner parameter.

    min_frac, max_frac: float, optional
        The minimum/maximum fraction of the axis width for the scale bar to
        extend. A value of 1 would allow the scale bar to extend across the
        entire axis width.  Only used for automatically calculating
        best-fit coeff and unit when neither is specified, otherwise
        disregarded.

    coord_system : string, optional
        This string defines the coordinate system of the coordinates of pos
        Valid coordinates are:

        - "data": 3D dataset coordinates
        - "plot": 2D coordinates defined by the actual plot limits
        - "axis": MPL axis coordinates: (0,0) is lower left; (1,1) is upper right
        - "figure": MPL figure coordinates: (0,0) is lower left, (1,1) is upper right

    text_args : dictionary, optional
        A dictionary of parameters to used to update the font_properties
        for the text in this callback.  For any property not set, it will
        use the defaults of the plot.  Thus one can modify the text size with
        ``text_args={'size':24}``

    size_bar_args : dictionary, optional
        A dictionary of parameters to be passed to the Matplotlib
        AnchoredSizeBar initializer.
        Defaults: ``{'pad': 0.25, 'sep': 5, 'borderpad': 1, 'color': 'w'}``

    draw_inset_box : boolean, optional
        Whether or not an inset box should be included around the scale bar.

    inset_box_args : dictionary, optional
        A dictionary of keyword arguments to be passed to the matplotlib Patch
        object that represents the inset box.
        Defaults: ``{'facecolor': 'black', 'linewidth': 3,
        'edgecolor': 'white', 'alpha': 0.5, 'boxstyle': 'square'}``

    scale_text_format : string, optional
        This specifies the format of the scalebar value assuming "scale" is the
        numerical value and "unit" is units of the scale (e.g. 'cm', 'kpc', etc.)
        The scale can be specified to arbitrary precision according to printf
        formatting codes. The format string must only specify "scale" and "units".
        Example: "Length = {scale:.2f} {units}". Default: "{scale} {units}"

    Example
    -------

    >>> import yt
    >>> ds = yt.load("Enzo_64/DD0020/data0020")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_scale()
    """

    _type_name = "scale"
    _supported_geometries = ("cartesian", "spectral_cube", "force")

    def __init__(
        self,
        *,
        corner="lower_right",
        coeff=None,
        unit=None,
        pos=None,
        max_frac=0.16,
        min_frac=0.015,
        coord_system="axis",
        text_args=None,
        size_bar_args=None,
        draw_inset_box=False,
        inset_box_args=None,
        scale_text_format="{scale} {units}",
    ):

        def_size_bar_args = {"pad": 0.05, "sep": 5, "borderpad": 1, "color": "w"}

        def_inset_box_args = {
            "facecolor": "black",
            "linewidth": 3,
            "edgecolor": "white",
            "alpha": 0.5,
            "boxstyle": "square",
        }

        # Set position based on corner argument.
        self.corner = corner
        self.coeff = coeff
        self.unit = unit
        self.pos = pos
        self.max_frac = max_frac
        self.min_frac = min_frac
        self.coord_system = coord_system
        self.scale_text_format = scale_text_format
        if size_bar_args is None:
            self.size_bar_args = def_size_bar_args
        else:
            self.size_bar_args = size_bar_args
        if inset_box_args is None:
            self.inset_box_args = def_inset_box_args
        else:
            self.inset_box_args = inset_box_args
        self.draw_inset_box = draw_inset_box
        if text_args is None:
            text_args = {}
        self.text_args = text_args

    def __call__(self, plot):
        from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

        # Callback only works for plots with axis ratios of 1
        xsize = plot.xlim[1] - plot.xlim[0]

        # Setting pos overrides corner argument
        if self.pos is None:
            if self.corner == "upper_left":
                self.pos = (0.11, 0.952)
            elif self.corner == "upper_right":
                self.pos = (0.89, 0.952)
            elif self.corner == "lower_left":
                self.pos = (0.11, 0.052)
            elif self.corner == "lower_right":
                self.pos = (0.89, 0.052)
            elif self.corner is None:
                self.pos = (0.5, 0.5)
            else:
                raise ValueError(
                    "Argument 'corner' must be set to "
                    "'upper_left', 'upper_right', 'lower_left', "
                    "'lower_right', or None"
                )

        # When identifying a best fit distance unit, do not allow scale marker
        # to be greater than max_frac fraction of xaxis or under min_frac
        # fraction of xaxis
        max_scale = self.max_frac * xsize
        min_scale = self.min_frac * xsize

        # If no units are set, pick something sensible.
        if self.unit is None:
            # User has set the axes units and supplied a coefficient.
            if plot._axes_unit_names is not None and self.coeff is not None:
                self.unit = plot._axes_unit_names[0]
            # Nothing provided; identify a best fit distance unit.
            else:
                min_scale = plot.ds.get_smallest_appropriate_unit(
                    min_scale, return_quantity=True
                )
                max_scale = plot.ds.get_smallest_appropriate_unit(
                    max_scale, return_quantity=True
                )
                if self.coeff is None:
                    self.coeff = max_scale.v
                self.unit = max_scale.units
        elif self.coeff is None:
            self.coeff = 1
        self.scale = plot.ds.quan(self.coeff, self.unit)
        text = self.scale_text_format.format(scale=int(self.coeff), units=self.unit)
        image_scale = (
            plot.frb.convert_distance_x(self.scale) / plot.frb.convert_distance_x(xsize)
        ).v
        size_vertical = self.size_bar_args.pop("size_vertical", 0.005 * plot.aspect)
        fontproperties = self.size_bar_args.pop(
            "fontproperties", plot.font_properties.copy()
        )
        frameon = self.size_bar_args.pop("frameon", self.draw_inset_box)
        # FontProperties instances use set_<property>() setter functions
        for key, val in self.text_args.items():
            setter_func = "set_" + key
            try:
                getattr(fontproperties, setter_func)(val)
            except AttributeError as e:
                raise AttributeError(
                    "Cannot set text_args keyword "
                    "to include '%s' because MPL's fontproperties object does "
                    "not contain function '%s'." % (key, setter_func)
                ) from e

        # this "anchors" the size bar to a box centered on self.pos in axis
        # coordinates
        self.size_bar_args["bbox_to_anchor"] = self.pos
        self.size_bar_args["bbox_transform"] = plot._axes.transAxes

        bar = AnchoredSizeBar(
            plot._axes.transAxes,
            image_scale,
            text,
            10,
            size_vertical=size_vertical,
            fontproperties=fontproperties,
            frameon=frameon,
            **self.size_bar_args,
        )

        bar.patch.set(**self.inset_box_args)

        plot._axes.add_artist(bar)

        return plot


class RayCallback(PlotCallback):
    """
    Adds a line representing the projected path of a ray across the plot.
    The ray can be either a YTOrthoRay, YTRay, or a LightRay object.
    annotate_ray() will properly account for periodic rays across the volume.
    If arrow is set to True, uses the MPL.pyplot.arrow function, otherwise
    uses the MPL.pyplot.plot function to plot a normal line.  Adjust
    plot_args accordingly.

    Parameters
    ----------

    ray : YTOrthoRay, YTRay, or LightRay
        Ray is the object that we want to include.  We overplot the projected
        trajectory of the ray.  If the object is a trident.LightRay
        object, it will only plot the segment of the LightRay that intersects
        the dataset currently displayed.

    arrow : boolean, optional
        Whether or not to place an arrowhead on the front of the ray to denote
        direction
        Default: False

    plot_args : dictionary, optional
        A dictionary of any arbitrary parameters to be passed to the Matplotlib
        line object.  Defaults: {'color':'white', 'linewidth':2}.

    Examples
    --------

    >>> # Overplot a ray and an ortho_ray object on a projection
    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> oray = ds.ortho_ray(1, (0.3, 0.4))  # orthoray down the y axis
    >>> ray = ds.ray((0.1, 0.2, 0.3), (0.6, 0.7, 0.8))  # arbitrary ray
    >>> p = yt.ProjectionPlot(ds, "z", "density")
    >>> p.annotate_ray(oray)
    >>> p.annotate_ray(ray)
    >>> p.save()

    >>> # Overplot a LightRay object on a projection
    >>> import yt
    >>> from trident import LightRay
    >>> ds = yt.load("enzo_cosmology_plus/RD0004/RD0004")
    >>> lr = LightRay(
    ...     "enzo_cosmology_plus/AMRCosmology.enzo", "Enzo", 0.0, 0.1, time_data=False
    ... )
    >>> lray = lr.make_light_ray(seed=1)
    >>> p = yt.ProjectionPlot(ds, "z", "density")
    >>> p.annotate_ray(lr)
    >>> p.save()

    """

    _type_name = "ray"
    _supported_geometries = ("cartesian", "spectral_cube", "force")

    def __init__(self, ray, *, arrow=False, plot_args=None, **kwargs):
        def_plot_args = {"color": "white", "linewidth": 2}
        self.ray = ray
        self.arrow = arrow
        if plot_args is not None:
            issue_deprecation_warning(
                "`plot_args` is deprecated. "
                "You can now pass arbitrary keyword arguments instead of a dictionary.",
                since="4.1.0",
            )
            plot_args = {**def_plot_args, **plot_args, **kwargs}
        else:
            plot_args = {**def_plot_args, **kwargs}
        self.plot_args = plot_args

    def _process_ray(self):
        """
        Get the start_coord and end_coord of a ray object
        """
        return (self.ray.start_point, self.ray.end_point)

    def _process_ortho_ray(self):
        """
        Get the start_coord and end_coord of an ortho_ray object
        """
        start_coord = self.ray.ds.domain_left_edge.copy()
        end_coord = self.ray.ds.domain_right_edge.copy()

        xax = self.ray.ds.coordinates.x_axis[self.ray.axis]
        yax = self.ray.ds.coordinates.y_axis[self.ray.axis]
        start_coord[xax] = end_coord[xax] = self.ray.coords[0]
        start_coord[yax] = end_coord[yax] = self.ray.coords[1]
        return (start_coord, end_coord)

    def _process_light_ray(self, plot):
        """
        Get the start_coord and end_coord of a LightRay object.
        Identify which of the sections of the LightRay is in the
        dataset that is currently being plotted.  If there is one, return the
        start and end of the corresponding ray segment
        """

        for ray_ds in self.ray.light_ray_solution:
            if ray_ds["unique_identifier"] == str(plot.ds.unique_identifier):
                start_coord = plot.ds.arr(ray_ds["start"])
                end_coord = plot.ds.arr(ray_ds["end"])
                return (start_coord, end_coord)
        # if no intersection between the plotted dataset and the LightRay
        # return a false tuple to pass to start_coord
        return ((False, False), (False, False))

    def __call__(self, plot):
        type_name = getattr(self.ray, "_type_name", None)

        if type_name == "ray":
            start_coord, end_coord = self._process_ray()

        elif type_name == "ortho_ray":
            start_coord, end_coord = self._process_ortho_ray()

        elif hasattr(self.ray, "light_ray_solution"):
            start_coord, end_coord = self._process_light_ray(plot)

        else:
            raise ValueError("ray must be a YTRay, YTOrthoRay, or LightRay object.")

        # if start_coord and end_coord are all False, it means no intersecting
        # ray segment with this plot.
        if not all(start_coord) and not all(end_coord):
            return plot

        # if possible, break periodic ray into non-periodic
        # segments and add each of them individually
        if any(plot.ds.periodicity):
            segments = periodic_ray(
                start_coord.to("code_length"),
                end_coord.to("code_length"),
                left=plot.ds.domain_left_edge.to("code_length"),
                right=plot.ds.domain_right_edge.to("code_length"),
            )
        else:
            segments = [[start_coord, end_coord]]

        # To assure that the last ray segment has an arrow if so desired
        # and all other ray segments are lines
        for segment in segments[:-1]:
            cb = LinePlotCallback(
                segment[0], segment[1], coord_system="data", **self.plot_args
            )
            cb(plot)
        segment = segments[-1]
        if self.arrow:
            cb = ArrowCallback(
                segment[1],
                starting_pos=segment[0],
                coord_system="data",
                **self.plot_args,
            )
        else:
            cb = LinePlotCallback(
                segment[0], segment[1], coord_system="data", **self.plot_args
            )
        cb(plot)
        return plot


class LineIntegralConvolutionCallback(PlotCallback):
    """
    Add the line integral convolution to the plot for vector fields
    visualization. Two component of vector fields needed to be provided
    (i.e., velocity_x and velocity_y, magnetic_field_x and magnetic_field_y).

    Parameters
    ----------

    field_x, field_y : string
        The names of two components of vector field which will be visualized

    texture : 2-d array with the same shape of image, optional
        Texture will be convolved when computing line integral convolution.
        A white noise background will be used as default.

    kernellen : float, optional
        The lens of kernel for convolution, which is the length over which the
        convolution will be performed. For longer kernellen, longer streamline
        structure will appear.

    lim : 2-element tuple, list, or array, optional
        The value of line integral convolution will be clipped to the range
        of lim, which applies upper and lower bounds to the values of line
        integral convolution and enhance the visibility of plots. Each element
        should be in the range of [0,1].

    cmap : string, optional
        The name of colormap for line integral convolution plot.

    alpha : float, optional
        The alpha value for line integral convolution plot.

    const_alpha : boolean, optional
        If set to False (by default), alpha will be weighted spatially by
        the values of line integral convolution; otherwise a constant value
        of the given alpha is used.

    Example
    -------

    >>> import yt
    >>> ds = yt.load("Enzo_64/DD0020/data0020")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_line_integral_convolution(
    ...     "velocity_x", "velocity_y", lim=(0.5, 0.65)
    ... )
    """

    _type_name = "line_integral_convolution"
    _supported_geometries = ("cartesian", "spectral_cube", "polar", "cylindrical")
    _incompatible_plot_types = ("LineIntegralConvolutionCallback",)

    def __init__(
        self,
        field_x,
        field_y,
        texture=None,
        kernellen=50.0,
        lim=(0.5, 0.6),
        cmap="binary",
        alpha=0.8,
        const_alpha=False,
    ):
        self.field_x = field_x
        self.field_y = field_y
        self.texture = texture
        self.kernellen = kernellen
        self.lim = lim
        self.cmap = cmap
        self.alpha = alpha
        self.const_alpha = const_alpha

    def __call__(self, plot):
        from matplotlib import cm

        bounds = self._physical_bounds(plot)
        extent = self._plot_bounds(plot)

        # We are feeding this size into the pixelizer, where it will properly
        # set it in reverse order
        nx = plot.raw_image_shape[1]
        ny = plot.raw_image_shape[0]
        pixX = plot.data.ds.coordinates.pixelize(
            plot.data.axis, plot.data, self.field_x, bounds, (nx, ny)
        )
        pixY = plot.data.ds.coordinates.pixelize(
            plot.data.axis, plot.data, self.field_y, bounds, (nx, ny)
        )

        vectors = np.concatenate((pixX[..., np.newaxis], pixY[..., np.newaxis]), axis=2)

        if self.texture is None:
            prng = np.random.RandomState(0x4D3D3D3)
            self.texture = prng.random_sample((nx, ny))
        elif self.texture.shape != (nx, ny):
            raise ValueError(
                "'texture' must have the same shape "
                "with that of output image (%d, %d)" % (nx, ny)
            )

        kernel = np.sin(np.arange(self.kernellen) * np.pi / self.kernellen)
        kernel = kernel.astype(np.double)

        lic_data = line_integral_convolution_2d(vectors, self.texture, kernel)
        lic_data = lic_data / lic_data.max()
        lic_data_clip = np.clip(lic_data, self.lim[0], self.lim[1])

        mask = ~(np.isfinite(pixX) & np.isfinite(pixY))
        lic_data[mask] = np.nan
        lic_data_clip[mask] = np.nan

        if plot._swap_axes:
            lic_data_clip = lic_data_clip.transpose()
            extent = (extent[2], extent[3], extent[0], extent[1])

        if self.const_alpha:
            plot._axes.imshow(
                lic_data_clip,
                extent=extent,
                cmap=self.cmap,
                alpha=self.alpha,
                origin="lower",
                aspect="auto",
            )
        else:
            lic_data_rgba = cm.ScalarMappable(norm=None, cmap=self.cmap).to_rgba(
                lic_data_clip
            )
            lic_data_clip_rescale = (lic_data_clip - self.lim[0]) / (
                self.lim[1] - self.lim[0]
            )
            lic_data_rgba[..., 3] = lic_data_clip_rescale * self.alpha
            plot._axes.imshow(
                lic_data_rgba,
                extent=extent,
                cmap=self.cmap,
                origin="lower",
                aspect="auto",
            )

        return plot


class CellEdgesCallback(PlotCallback):
    """
    Annotate cell edges.  This is done through a second call to pixelize, where
    the distance from a pixel to a cell boundary in pixels is compared against
    the `line_width` argument.  The secondary image is colored as `color` and
    overlaid with the `alpha` value.

    Parameters
    ----------
    line_width : float
        The width of the cell edge lines in normalized units relative to the
        size of the longest axis.  Default is 1% of the size of the smallest
        axis.
    alpha : float
        When the second image is overlaid, it will have this level of alpha
        transparency.  Default is 1.0 (fully-opaque).
    color : tuple of three floats or matplotlib color name
        This is the color of the cell edge values.  It defaults to black.

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> s = yt.SlicePlot(ds, "z", "density")
    >>> s.annotate_cell_edges()
    >>> s.save()
    """

    _type_name = "cell_edges"
    _supported_geometries = ("cartesian", "spectral_cube", "cylindrical")

    def __init__(self, line_width=0.002, alpha=1.0, color="black"):
        from matplotlib.colors import ColorConverter

        conv = ColorConverter()
        self.line_width = line_width
        self.alpha = alpha
        self.color = (np.array(conv.to_rgb(color)) * 255).astype("uint8")

    def __call__(self, plot):
        if plot.data.ds.geometry == "cylindrical" and plot.data.ds.dimensionality == 3:
            raise NotImplementedError(
                "Cell edge annotation is only supported for \
                for 2D cylindrical geometry, not 3D"
            )
        x0, x1, y0, y1 = self._physical_bounds(plot)
        nx = plot.raw_image_shape[1]
        ny = plot.raw_image_shape[0]
        aspect = float((y1 - y0) / (x1 - x0))
        pixel_aspect = float(ny) / nx
        relative_aspect = pixel_aspect / aspect
        if relative_aspect > 1:
            nx = int(nx / relative_aspect)
        else:
            ny = int(ny * relative_aspect)
        if aspect > 1:
            if nx < 1600:
                nx = int(1600.0 / nx * ny)
                ny = 1600
            long_axis = ny
        else:
            if ny < 1600:
                nx = int(1600.0 / ny * nx)
                ny = 1600
            long_axis = nx
        line_width = max(self.line_width * long_axis, 1.0)
        im = np.zeros((ny, nx), dtype="f8")
        pixelize_cartesian(
            im,
            plot.data["px"],
            plot.data["py"],
            plot.data["pdx"],
            plot.data["pdy"],
            plot.data["px"],  # dummy field
            (x0, x1, y0, y1),
            line_width=line_width,
        )
        # New image:
        im_buffer = np.zeros((ny, nx, 4), dtype="uint8")
        im_buffer[im > 0, 3] = 255
        im_buffer[im > 0, :3] = self.color

        extent = self._plot_bounds(plot)
        if plot._swap_axes:
            im_buffer = im_buffer.transpose((1, 0, 2))
        plot._axes.imshow(
            im_buffer,
            origin="lower",
            interpolation="bilinear",
            extent=extent,
            alpha=self.alpha,
        )
        self._set_plot_limits(plot, extent)
