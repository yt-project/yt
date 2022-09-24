import abc
from collections import defaultdict
from numbers import Number
from typing import List, Optional, Type, Union

import matplotlib
import numpy as np
from matplotlib.colors import Normalize
from more_itertools import always_iterable
from unyt.exceptions import UnitConversionError

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.data_objects.image_array import ImageArray
from yt.frontends.ytdata.data_structures import YTSpatialPlotDataset
from yt.funcs import (
    fix_axis,
    fix_unitary,
    is_sequence,
    iter_fields,
    mylog,
    obj_length,
    validate_moment,
)
from yt.units.unit_object import Unit  # type: ignore
from yt.units.unit_registry import UnitParseError  # type: ignore
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.exceptions import (
    YTCannotParseUnitDisplayName,
    YTDataTypeUnsupported,
    YTInvalidFieldType,
    YTPlotCallbackError,
    YTUnitNotRecognized,
)
from yt.utilities.math_utils import ortho_find
from yt.utilities.orientation import Orientation
from yt.visualization._handlers import ColorbarHandler, NormHandler
from yt.visualization.base_plot_types import CallbackWrapper, ImagePlotMPL

from ._commons import _swap_axes_extents, get_default_from_config
from .fixed_resolution import (
    FixedResolutionBuffer,
    OffAxisProjectionFixedResolutionBuffer,
)
from .geo_plot_utils import get_mpl_transform
from .plot_container import (
    ImagePlotContainer,
    invalidate_data,
    invalidate_figure,
    invalidate_plot,
)

import sys  # isort: skip

if sys.version_info >= (3, 10):
    pass
else:
    from yt._maintenance.backports import zip


def get_window_parameters(axis, center, width, ds):
    width = ds.coordinates.sanitize_width(axis, width, None)
    center, display_center = ds.coordinates.sanitize_center(center, axis)
    xax = ds.coordinates.x_axis[axis]
    yax = ds.coordinates.y_axis[axis]
    bounds = (
        display_center[xax] - width[0] / 2,
        display_center[xax] + width[0] / 2,
        display_center[yax] - width[1] / 2,
        display_center[yax] + width[1] / 2,
    )
    return (bounds, center, display_center)


def get_oblique_window_parameters(normal, center, width, ds, depth=None):
    display_center, center = ds.coordinates.sanitize_center(center, 4)
    width = ds.coordinates.sanitize_width(normal, width, depth)

    if len(width) == 2:
        # Transforming to the cutting plane coordinate system
        center = (center - ds.domain_left_edge) / ds.domain_width - 0.5
        (normal, perp1, perp2) = ortho_find(normal)
        mat = np.transpose(np.column_stack((perp1, perp2, normal)))
        center = np.dot(mat, center)

    w = tuple(el.in_units("code_length") for el in width)
    bounds = tuple(((2 * (i % 2)) - 1) * w[i // 2] / 2 for i in range(len(w) * 2))

    return (bounds, center)


def get_axes_unit(width, ds):
    r"""
    Infers the axes unit names from the input width specification
    """
    if ds.no_cgs_equiv_length:
        return ("code_length",) * 2
    if is_sequence(width):
        if isinstance(width[1], str):
            axes_unit = (width[1], width[1])
        elif is_sequence(width[1]):
            axes_unit = (width[0][1], width[1][1])
        elif isinstance(width[0], YTArray):
            axes_unit = (str(width[0].units), str(width[1].units))
        else:
            axes_unit = None
    else:
        if isinstance(width, YTArray):
            axes_unit = (str(width.units), str(width.units))
        else:
            axes_unit = None
    return axes_unit


def validate_mesh_fields(data_source, fields):
    # this check doesn't make sense for ytdata plot datasets, which
    # load mesh data as a particle field but nonetheless can still
    # make plots with it
    if isinstance(data_source.ds, YTSpatialPlotDataset):
        return
    canonical_fields = data_source._determine_fields(fields)
    invalid_fields = []
    for field in canonical_fields:
        finfo = data_source.ds.field_info[field]
        if finfo.sampling_type == "particle":
            if not hasattr(data_source.ds, "_sph_ptypes"):
                pass
            elif finfo.is_sph_field:
                continue
            invalid_fields.append(field)

    if len(invalid_fields) > 0:
        raise YTInvalidFieldType(invalid_fields)


class PlotWindow(ImagePlotContainer, abc.ABC):
    r"""
    A plotting mechanism based around the concept of a window into a
    data source. It can have arbitrary fields, each of which will be
    centered on the same viewpoint, but will have individual zlimits.

    The data and plot are updated separately, and each can be
    invalidated as the object is modified.

    Data is handled by a FixedResolutionBuffer object.

    Parameters
    ----------

    data_source :
        :class:`yt.data_objects.selection_objects.base_objects.YTSelectionContainer2D`
        This is the source to be pixelized, which can be a projection,
        slice, or a cutting plane.
    bounds : sequence of floats
        Bounds are the min and max in the image plane that we want our
        image to cover.  It's in the order of (xmin, xmax, ymin, ymax),
        where the coordinates are all in the appropriate code units.
    buff_size : sequence of ints
        The size of the image to generate.
    antialias : boolean
        This can be true or false.  It determines whether or not sub-pixel
        rendering is used during data deposition.
    window_size : float
        The size of the window on the longest axis (in units of inches),
        including the margins but not the colorbar.
    right_handed : boolean
        Depreceated, please use flip_horizontal callback.
        Whether the implicit east vector for the image generated is set to make a right
        handed coordinate system with a north vector and the normal vector, the
        direction of the 'window' into the data.

    """

    def __init__(
        self,
        data_source,
        bounds,
        buff_size=(800, 800),
        antialias=True,
        periodic=True,
        origin="center-window",
        oblique=False,
        right_handed=True,
        window_size=8.0,
        fields=None,
        fontsize=18,
        aspect=None,
        setup=False,
        *,
        geometry="cartesian",
    ):

        # axis manipulation operations are callback-only:
        self._swap_axes_input = False
        self._flip_vertical = False
        self._flip_horizontal = False

        self.center = None
        self._periodic = periodic
        self.oblique = oblique
        self._right_handed = _check_right_handed(right_handed, self._flip_horizontal)
        self._equivalencies = defaultdict(lambda: (None, {}))
        self.buff_size = buff_size
        self.antialias = antialias
        self._axes_unit_names = None
        self._transform = None
        self._projection = None

        self.aspect = aspect
        skip = list(FixedResolutionBuffer._exclude_fields) + data_source._key_fields

        fields = list(iter_fields(fields))
        self.override_fields = list(set(fields).intersection(set(skip)))
        self.fields = [f for f in fields if f not in skip]
        super().__init__(data_source, window_size, fontsize)

        self._set_window(bounds)  # this automatically updates the data and plot
        if (
            geometry
            in (
                "spherical",
                "cylindrical",
                "geographic",
                "internal_geographic",
                "polar",
            )
            and origin != "native"
        ):
            mylog.info("Setting origin='native' for %s geometry.", geometry)
            origin = "native"

        self.origin = origin
        if self.data_source.center is not None and not oblique:
            ax = self.data_source.axis
            xax = self.ds.coordinates.x_axis[ax]
            yax = self.ds.coordinates.y_axis[ax]
            center, display_center = self.ds.coordinates.sanitize_center(
                self.data_source.center, ax
            )
            center = [display_center[xax], display_center[yax]]
            self.set_center(center)

            axname = self.ds.coordinates.axis_name[ax]
            transform = self.ds.coordinates.data_transform[axname]
            projection = self.ds.coordinates.data_projection[axname]
            self._projection = get_mpl_transform(projection)
            self._transform = get_mpl_transform(transform)

        self._setup_plots()

        for field in self.data_source._determine_fields(self.fields):
            finfo = self.data_source.ds._get_field_info(*field)
            pnh = self.plots[field].norm_handler
            if finfo.take_log is False:
                # take_log can be `None` so we explicitly compare against a boolean
                pnh.norm_type = Normalize
            else:
                # do nothing, the norm handler is responsible for
                # determining a viable norm, and defaults to LogNorm/SymLogNorm
                pass

            # override from user configuration if any
            log, linthresh = get_default_from_config(
                self.data_source,
                field=field,
                keys=["log", "linthresh"],
                defaults=[None, None],
            )
            if linthresh is not None:
                self.set_log(field, linthresh=linthresh)
            elif log is not None:
                self.set_log(field, log)

    def __iter__(self):
        for ds in self.ts:
            mylog.warning("Switching to %s", ds)
            self._switch_ds(ds)
            yield self

    def piter(self, *args, **kwargs):
        for ds in self.ts.piter(*args, **kwargs):
            self._switch_ds(ds)
            yield self

    _frb = None

    @property
    def frb(self):
        # Force the regeneration of the fixed resolution buffer
        # * if there's none
        # * if the data has been invalidated
        # * if the frb has been inalidated
        if not self._data_valid:
            self._recreate_frb()
        return self._frb

    @frb.setter
    def frb(self, value):
        self._frb = value
        self._data_valid = True

    @frb.deleter
    def frb(self):
        del self._frb
        self._frb = None

    def _recreate_frb(self):
        old_fields = None
        old_filters = []
        # If we are regenerating an frb, we want to know what fields we had before
        if self._frb is not None:
            old_fields = list(self._frb.data.keys())
            old_units = [_.units for _ in self._frb.data.values()]
            old_filters = self._frb._filters
        # Set the bounds
        if hasattr(self, "zlim"):
            # Support OffAxisProjectionPlot and OffAxisSlicePlot
            bounds = self.xlim + self.ylim + self.zlim
        else:
            bounds = self.xlim + self.ylim

        # Generate the FRB
        self.frb = self._frb_generator(
            self.data_source,
            bounds,
            self.buff_size,
            self.antialias,
            periodic=self._periodic,
            filters=old_filters,
        )

        # At this point the frb has the valid bounds, size, aliasing, etc.
        if old_fields is not None:
            # Restore the old fields
            for key, units in zip(old_fields, old_units):
                self._frb.render(key)
                equiv = self._equivalencies[key]
                if equiv[0] is None:
                    self._frb[key].convert_to_units(units)
                else:
                    self.frb.set_unit(key, units, equiv[0], equiv[1])

        # Restore the override fields
        for key in self.override_fields:
            self._frb.render(key)

    @property
    def _has_swapped_axes(self):
        # note: we always run the validations here in case the states of
        # the conflicting attributes have changed.
        return self._validate_swap_axes(self._swap_axes_input)

    @invalidate_data
    def swap_axes(self):
        # toggles the swap_axes behavior
        new_swap_value = not self._swap_axes_input
        # note: we also validate here to catch invalid states immediately, even
        # though we validate on accessing the attribute in `_has_swapped_axes`.
        self._swap_axes_input = self._validate_swap_axes(new_swap_value)
        return self

    def _validate_swap_axes(self, swap_value: bool) -> bool:
        if swap_value and (self._transform or self._projection):
            mylog.warning("Cannot swap axes due to transform or projection")
            return False
        return swap_value

    @property
    def width(self):
        Wx = self.xlim[1] - self.xlim[0]
        Wy = self.ylim[1] - self.ylim[0]
        return (Wx, Wy)

    @property
    def bounds(self):
        return self.xlim + self.ylim

    @invalidate_data
    def zoom(self, factor):
        r"""This zooms the window by *factor* > 0.
        - zoom out with *factor* < 1
        - zoom in with *factor* > 1

        Parameters
        ----------
        factor : float
            multiplier for the current width

        """
        if factor <= 0:
            raise ValueError("Only positive zooming factors are meaningful.")
        Wx, Wy = self.width
        centerx = self.xlim[0] + Wx * 0.5
        centery = self.ylim[0] + Wy * 0.5
        nWx, nWy = Wx / factor, Wy / factor
        self.xlim = (centerx - nWx * 0.5, centerx + nWx * 0.5)
        self.ylim = (centery - nWy * 0.5, centery + nWy * 0.5)
        return self

    @invalidate_data
    def pan(self, deltas):
        r"""Pan the image by specifying absolute code unit coordinate deltas.

        Parameters
        ----------
        deltas : Two-element sequence of floats, quantities, or (float, unit)
                 tuples.

            (delta_x, delta_y).  If a unit is not supplied the unit is assumed
            to be code_length.

        """
        if len(deltas) != 2:
            raise TypeError(
                f"The pan function accepts a two-element sequence.\nReceived {deltas}."
            )
        if isinstance(deltas[0], Number) and isinstance(deltas[1], Number):
            deltas = (
                self.ds.quan(deltas[0], "code_length"),
                self.ds.quan(deltas[1], "code_length"),
            )
        elif isinstance(deltas[0], tuple) and isinstance(deltas[1], tuple):
            deltas = (
                self.ds.quan(deltas[0][0], deltas[0][1]),
                self.ds.quan(deltas[1][0], deltas[1][1]),
            )
        elif isinstance(deltas[0], YTQuantity) and isinstance(deltas[1], YTQuantity):
            pass
        else:
            raise TypeError(
                "The arguments of the pan function must be a sequence of floats,\n"
                "quantities, or (float, unit) tuples. Received %s." % (deltas,)
            )
        self.xlim = (self.xlim[0] + deltas[0], self.xlim[1] + deltas[0])
        self.ylim = (self.ylim[0] + deltas[1], self.ylim[1] + deltas[1])
        return self

    @invalidate_data
    def pan_rel(self, deltas):
        r"""Pan the image by specifying relative deltas, to the FOV.

        Parameters
        ----------
        deltas : sequence of floats
            (delta_x, delta_y) in *relative* code unit coordinates

        """
        Wx, Wy = self.width
        self.xlim = (self.xlim[0] + Wx * deltas[0], self.xlim[1] + Wx * deltas[0])
        self.ylim = (self.ylim[0] + Wy * deltas[1], self.ylim[1] + Wy * deltas[1])
        return self

    @invalidate_plot
    def set_unit(self, field, new_unit, equivalency=None, equivalency_kwargs=None):
        """Sets a new unit for the requested field

        parameters
        ----------
        field : string or field tuple
           The name of the field that is to be changed.

        new_unit : string or Unit object

        equivalency : string, optional
           If set, the equivalency to use to convert the current units to
           the new requested unit. If None, the unit conversion will be done
           without an equivalency

        equivalency_kwargs : string, optional
           Keyword arguments to be passed to the equivalency. Only used if
           ``equivalency`` is set.
        """
        for f, u in zip(iter_fields(field), always_iterable(new_unit), strict=True):
            self.frb.set_unit(f, u, equivalency, equivalency_kwargs)
            self._equivalencies[f] = (equivalency, equivalency_kwargs)
            pnh = self.plots[f].norm_handler
            pnh.display_units = u
        return self

    @invalidate_plot
    def set_origin(self, origin):
        """Set the plot origin.

        Parameters
        ----------
        origin : string or length 1, 2, or 3 sequence.
           The location of the origin of the plot coordinate system. This
           is typically represented by a '-' separated string or a tuple of
           strings. In the first index the y-location is given by 'lower',
           'upper', or 'center'. The second index is the x-location, given as
           'left', 'right', or 'center'. Finally, whether the origin is
           applied in 'domain' space, plot 'window' space or 'native'
           simulation coordinate system is given. For example, both
           'upper-right-domain' and ['upper', 'right', 'domain'] place the
           origin in the upper right hand corner of domain space. If x or y
           are not given, a value is inferred. For instance, 'left-domain'
           corresponds to the lower-left hand corner of the simulation domain,
           'center-domain' corresponds to the center of the simulation domain,
           or 'center-window' for the center of the plot window. In the event
           that none of these options place the origin in a desired location,
           a sequence of tuples and a string specifying the
           coordinate space can be given. If plain numeric types are input,
           units of `code_length` are assumed. Further examples:

        ===============================================  ===============================
        format                                           example
        ===============================================  ===============================
        '{space}'                                        'domain'
        '{xloc}-{space}'                                 'left-window'
        '{yloc}-{space}'                                 'upper-domain'
        '{yloc}-{xloc}-{space}'                          'lower-right-window'
        ('{space}',)                                     ('window',)
        ('{xloc}', '{space}')                            ('right', 'domain')
        ('{yloc}', '{space}')                            ('lower', 'window')
        ('{yloc}', '{xloc}', '{space}')                  ('lower', 'right', 'window')
        ((yloc, '{unit}'), (xloc, '{unit}'), '{space}')  ((0, 'm'), (.4, 'm'), 'window')
        (xloc, yloc, '{space}')                          (0.23, 0.5, 'domain')
        ===============================================  ===============================
        """
        self.origin = origin
        return self

    @invalidate_plot
    @invalidate_figure
    def set_mpl_projection(self, mpl_proj):
        r"""
        Set the matplotlib projection type with a cartopy transform function

        Given a string or a tuple argument, this will project the data onto
        the plot axes with the chosen transform function.

        Assumes that the underlying data has a PlateCarree transform type.

        To annotate the plot with coastlines or other annotations,
        `_setup_plots()` will need to be called after this function
        to make the axes available for annotation.

        Parameters
        ----------

        mpl_proj : string or tuple
           if passed as a string, mpl_proj is the specified projection type,
           if passed as a tuple, then tuple will take the form of
           ``("ProjectionType", (args))`` or ``("ProjectionType", (args), {kwargs})``
           Valid projection type options include:
           'PlateCarree', 'LambertConformal', 'LabmbertCylindrical',
           'Mercator', 'Miller', 'Mollweide', 'Orthographic',
           'Robinson', 'Stereographic', 'TransverseMercator',
           'InterruptedGoodeHomolosine', 'RotatedPole', 'OGSB',
           'EuroPP', 'Geostationary', 'Gnomonic', 'NorthPolarStereo',
           'OSNI', 'SouthPolarStereo', 'AlbersEqualArea',
           'AzimuthalEquidistant', 'Sinusoidal', 'UTM',
           'NearsidePerspective', 'LambertAzimuthalEqualArea'

        Examples
        --------

        This will create a Mollweide projection using Mollweide default values
        and annotate it with coastlines

        >>> import yt
        >>> ds = yt.load("")
        >>> p = yt.SlicePlot(ds, "altitude", "AIRDENS")
        >>> p.set_mpl_projection("AIRDENS", "Mollweide")
        >>> p._setup_plots()
        >>> p.plots["AIRDENS"].axes.coastlines()
        >>> p.show()

        This will move the PlateCarree central longitude to 90 degrees and
        annotate with coastlines.

        >>> import yt
        >>> ds = yt.load("")
        >>> p = yt.SlicePlot(ds, "altitude", "AIRDENS")
        >>> p.set_mpl_projection(
        ...     "AIRDENS", ("PlateCarree", (), {"central_longitude": 90, "globe": None})
        ... )
        >>> p._setup_plots()
        >>> p.plots["AIRDENS"].axes.set_global()
        >>> p.plots["AIRDENS"].axes.coastlines()
        >>> p.show()


        This will create a RoatatedPole projection with the unrotated pole
        position at 37.5 degrees latitude and 177.5 degrees longitude by
        passing them in as args.


        >>> import yt
        >>> ds = yt.load("")
        >>> p = yt.SlicePlot(ds, "altitude", "AIRDENS")
        >>> p.set_mpl_projection("RotatedPole", (177.5, 37.5))
        >>> p._setup_plots()
        >>> p.plots["AIRDENS"].axes.set_global()
        >>> p.plots["AIRDENS"].axes.coastlines()
        >>> p.show()

        This will create a RoatatedPole projection with the unrotated pole
        position at 37.5 degrees latitude and 177.5 degrees longitude by
        passing them in as kwargs.

        >>> import yt
        >>> ds = yt.load("")
        >>> p = yt.SlicePlot(ds, "altitude", "AIRDENS")
        >>> p.set_mpl_projection(
        ...     ("RotatedPole", (), {"pole_latitude": 37.5, "pole_longitude": 177.5})
        ... )
        >>> p._setup_plots()
        >>> p.plots["AIRDENS"].axes.set_global()
        >>> p.plots["AIRDENS"].axes.coastlines()
        >>> p.show()

        """

        self._projection = get_mpl_transform(mpl_proj)
        axname = self.ds.coordinates.axis_name[self.data_source.axis]
        transform = self.ds.coordinates.data_transform[axname]
        self._transform = get_mpl_transform(transform)
        return self

    @invalidate_data
    def _set_window(self, bounds):
        """Set the bounds of the plot window.
        This is normally only called internally, see set_width.


        Parameters
        ----------

        bounds : a four element sequence of floats
            The x and y bounds, in the format (x0, x1, y0, y1)

        """
        if self.center is not None:
            dx = bounds[1] - bounds[0]
            dy = bounds[3] - bounds[2]
            self.xlim = (self.center[0] - dx / 2.0, self.center[0] + dx / 2.0)
            self.ylim = (self.center[1] - dy / 2.0, self.center[1] + dy / 2.0)
        else:
            self.xlim = tuple(bounds[0:2])
            self.ylim = tuple(bounds[2:4])
            if len(bounds) == 6:
                # Support OffAxisProjectionPlot and OffAxisSlicePlot
                self.zlim = tuple(bounds[4:6])
        mylog.info("xlim = %f %f", self.xlim[0], self.xlim[1])
        mylog.info("ylim = %f %f", self.ylim[0], self.ylim[1])
        if hasattr(self, "zlim"):
            mylog.info("zlim = %f %f", self.zlim[0], self.zlim[1])

    @invalidate_data
    def set_width(self, width, unit=None):
        """set the width of the plot window

        parameters
        ----------
        width : float, array of floats, (float, unit) tuple, or tuple of
                (float, unit) tuples.

             Width can have four different formats to support windows with
             variable x and y widths.  They are:

             ==================================     =======================
             format                                 example
             ==================================     =======================
             (float, string)                        (10,'kpc')
             ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
             float                                  0.2
             (float, float)                         (0.2, 0.3)
             ==================================     =======================

             For example, (10, 'kpc') requests a plot window that is 10
             kiloparsecs wide in the x and y directions,
             ((10,'kpc'),(15,'kpc')) requests a window that is 10 kiloparsecs
             wide along the x axis and 15 kiloparsecs wide along the y axis.
             In the other two examples, code units are assumed, for example
             (0.2, 0.3) requests a plot that has an x width of 0.2 and a y
             width of 0.3 in code units.  If units are provided the resulting
             plot axis labels will use the supplied units.
        unit : str
             the unit the width has been specified in. If width is a tuple, this
             argument is ignored. Defaults to code units.
        """

        if isinstance(width, Number):
            if unit is None:
                width = (width, "code_length")
            else:
                width = (width, fix_unitary(unit))

        axes_unit = get_axes_unit(width, self.ds)

        width = self.ds.coordinates.sanitize_width(self.frb.axis, width, None)

        centerx = (self.xlim[1] + self.xlim[0]) / 2.0
        centery = (self.ylim[1] + self.ylim[0]) / 2.0

        self.xlim = (centerx - width[0] / 2, centerx + width[0] / 2)
        self.ylim = (centery - width[1] / 2, centery + width[1] / 2)

        if hasattr(self, "zlim"):
            centerz = (self.zlim[1] + self.zlim[0]) / 2.0
            mw = self.ds.arr(width).max()
            self.zlim = (centerz - mw / 2.0, centerz + mw / 2.0)

        self.set_axes_unit(axes_unit)

        return self

    @invalidate_data
    def set_center(self, new_center, unit="code_length"):
        """Sets a new center for the plot window

        parameters
        ----------
        new_center : two element sequence of floats
            The coordinates of the new center of the image in the
            coordinate system defined by the plot axes. If the unit
            keyword is not specified, the coordinates are assumed to
            be in code units.

        unit : string
            The name of the unit new_center is given in.  If new_center is a
            YTArray or tuple of YTQuantities, this keyword is ignored.

        """
        error = RuntimeError(
            "\n"
            "new_center must be a two-element list or tuple of floats \n"
            "corresponding to a coordinate in the plot relative to \n"
            "the plot coordinate system.\n"
        )
        if new_center is None:
            self.center = None
        elif is_sequence(new_center):
            if len(new_center) != 2:
                raise error
            for el in new_center:
                if not isinstance(el, Number) and not isinstance(el, YTQuantity):
                    raise error
            if isinstance(new_center[0], Number):
                new_center = [self.ds.quan(c, unit) for c in new_center]
            self.center = new_center
        else:
            raise error
        self._set_window(self.bounds)
        return self

    @invalidate_data
    def set_antialias(self, aa):
        """Turn antialiasing on or off.

        parameters
        ----------
        aa : boolean
        """
        self.antialias = aa

    @invalidate_data
    def set_buff_size(self, size):
        """Sets a new buffer size for the fixed resolution buffer

        parameters
        ----------
        size : int or two element sequence of ints
            The number of data elements in the buffer on the x and y axes.
            If a scalar is provided,  then the buffer is assumed to be square.
        """
        if is_sequence(size):
            self.buff_size = size
        else:
            self.buff_size = (size, size)
        return self

    @invalidate_plot
    def set_axes_unit(self, unit_name):
        r"""Set the unit for display on the x and y axes of the image.

        Parameters
        ----------
        unit_name : string or two element tuple of strings
            A unit, available for conversion in the dataset, that the
            image extents will be displayed in.  If set to None, any previous
            units will be reset.  If the unit is None, the default is chosen.
            If unit_name is '1', 'u', or 'unitary', it will not display the
            units, and only show the axes name. If unit_name is a tuple, the
            first element is assumed to be the unit for the x axis and the
            second element the unit for the y axis.

        Raises
        ------
        YTUnitNotRecognized
            If the unit is not known, this will be raised.

        Examples
        --------

        >>> from yt import load
        >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> p = ProjectionPlot(ds, "y", "Density")
        >>> p.set_axes_unit("kpc")

        """
        # blind except because it could be in conversion_factors or units
        if unit_name is not None:
            if isinstance(unit_name, str):
                unit_name = (unit_name, unit_name)
            for un in unit_name:
                try:
                    self.ds.length_unit.in_units(un)
                except (UnitConversionError, UnitParseError) as e:
                    raise YTUnitNotRecognized(un) from e
        self._axes_unit_names = unit_name
        return self

    @invalidate_plot
    def toggle_right_handed(self):
        issue_deprecation_warning(
            "the toggle_right_handed method is deprecated, use `.flip_horizontal()` instead.",
            since="4.1.0",
            removal="4.2.0",
        )
        self.flip_horizontal()
        return self

    @invalidate_plot
    def flip_horizontal(self):
        """
        inverts the horizontal axis (the image's abscissa)
        """
        self._flip_horizontal = not self._flip_horizontal
        self._right_handed = (
            not self._right_handed
        )  # keep in sync until full depreciation
        return self

    @invalidate_plot
    def flip_vertical(self):
        """
        inverts the vertical axis (the image's ordinate)
        """
        self._flip_vertical = not self._flip_vertical
        return self

    def to_fits_data(self, fields=None, other_keys=None, length_unit=None, **kwargs):
        r"""Export the fields in this PlotWindow instance
        to a FITSImageData instance.

        This will export a set of FITS images of either the fields specified
        or all the fields already in the object.

        Parameters
        ----------
        fields : list of strings
            These fields will be pixelized and output. If "None", the keys of
            the FRB will be used.
        other_keys : dictionary, optional
            A set of header keys and values to write into the FITS header.
        length_unit : string, optional
            the length units that the coordinates are written in. The default
            is to use the default length unit of the dataset.
        """
        return self.frb.to_fits_data(
            fields=fields, other_keys=other_keys, length_unit=length_unit, **kwargs
        )


class PWViewerMPL(PlotWindow):
    """Viewer using matplotlib as a backend via the WindowPlotMPL."""

    _current_field = None
    _frb_generator: Optional[Type[FixedResolutionBuffer]] = None
    _plot_type: Optional[str] = None

    def __init__(self, *args, **kwargs):
        if self._frb_generator is None:
            self._frb_generator = kwargs.pop("frb_generator")
        if self._plot_type is None:
            self._plot_type = kwargs.pop("plot_type")
        self._splat_color = kwargs.pop("splat_color", None)
        self._frb: Optional[FixedResolutionBuffer] = None
        PlotWindow.__init__(self, *args, **kwargs)

        # import type here to avoid import cycles
        # note that this import statement is actually crucial at runtime:
        # the filter methods for the present class are defined only when
        # fixed_resolution_filters is imported, so we need to guarantee
        # that it happens no later than instanciation
        from yt.visualization.plot_modifications import PlotCallback

        self._callbacks: List[PlotCallback] = []

    @property
    def _data_valid(self) -> bool:
        return self._frb is not None and self._frb._data_valid

    @_data_valid.setter
    def _data_valid(self, value):
        if self._frb is None:
            # we delegate the (in)validation responsability to the FRB
            # if we don't have one yet, we can exit without doing anything
            return
        else:
            self._frb._data_valid = value

    def _setup_origin(self):
        origin = self.origin
        axis_index = self.data_source.axis
        xc = None
        yc = None

        if isinstance(origin, str):
            origin = tuple(origin.split("-"))

        if len(origin) > 3:
            raise ValueError(
                "Invalid origin argument with too many elements; "
                f"expected 1, 2 or 3 elements, got {self.origin!r}, counting {len(origin)} elements. "
                "Use '-' as a separator for string arguments."
            )

        if len(origin) == 1:
            coord_system = origin[0]
            if coord_system not in ("window", "domain", "native"):
                raise ValueError(
                    "Invalid origin argument. "
                    "Single element specification must be 'window', 'domain', or 'native'. "
                    f"Got {self.origin!r}"
                )
            origin = ("lower", "left", coord_system)

        elif len(origin) == 2:
            err_msg = "Invalid origin argument. Using 2 elements:\n"

            if origin[0] in ("left", "right", "center"):
                o0map = {"left": "lower", "right": "upper", "center": "center"}
                origin = (o0map[origin[0]],) + origin
            elif origin[0] in ("lower", "upper"):
                origin = (origin[0], "center", origin[-1])
            else:
                err_msg += " - the first one must be 'left', 'right', 'lower', 'upper' or 'center'\n"

            if origin[-1] not in ("window", "domain", "native"):
                err_msg += " - the second one must be 'window', 'domain', or 'native'\n"

            if len(err_msg.split("\n")) > 2:
                err_msg += f"Got {self.origin!r}"
                raise ValueError(err_msg)

        elif len(origin) == 3:
            err_msg = "Invalid origin argument. Using 3 elements:\n"
            if isinstance(origin[0], (int, float)):
                xc = self.ds.quan(origin[0], "code_length")
            elif isinstance(origin[0], tuple):
                xc = self.ds.quan(*origin[0])
            elif origin[0] not in ("lower", "upper", "center"):
                err_msg += " - the first one must be 'lower', 'upper' or 'center' or a distance\n"

            if isinstance(origin[1], (int, float)):
                yc = self.ds.quan(origin[1], "code_length")
            elif isinstance(origin[1], tuple):
                yc = self.ds.quan(*origin[1])
            elif origin[1] not in ("left", "right", "center"):
                err_msg += " - the second one must be 'left', 'right', 'center' or a distance\n"

            if origin[-1] not in ("window", "domain", "native"):
                err_msg += " - the third one must be 'window', 'domain', or 'native'\n"

            if len(err_msg.split("\n")) > 2:
                err_msg += f"Got {self.origin!r}"
                raise ValueError(err_msg)

        assert not isinstance(origin, str)
        assert len(origin) == 3
        assert origin[2] in ("window", "domain", "native")

        if origin[2] == "window":
            xllim, xrlim = self.xlim
            yllim, yrlim = self.ylim
        elif origin[2] == "domain":
            xax = self.ds.coordinates.x_axis[axis_index]
            yax = self.ds.coordinates.y_axis[axis_index]
            xllim = self.ds.domain_left_edge[xax]
            xrlim = self.ds.domain_right_edge[xax]
            yllim = self.ds.domain_left_edge[yax]
            yrlim = self.ds.domain_right_edge[yax]
        elif origin[2] == "native":
            return (self.ds.quan(0.0, "code_length"), self.ds.quan(0.0, "code_length"))

        if xc is None and yc is None:
            assert origin[0] in ("lower", "upper", "center")
            assert origin[1] in ("left", "right", "center")

            if origin[0] == "lower":
                yc = yllim
            elif origin[0] == "upper":
                yc = yrlim
            elif origin[0] == "center":
                yc = (yllim + yrlim) / 2.0

            if origin[1] == "left":
                xc = xllim
            elif origin[1] == "right":
                xc = xrlim
            elif origin[1] == "center":
                xc = (xllim + xrlim) / 2.0

        x_in_bounds = xc >= xllim and xc <= xrlim
        y_in_bounds = yc >= yllim and yc <= yrlim

        if not x_in_bounds and not y_in_bounds:
            raise ValueError(
                "origin inputs not in bounds of specified coordinate system domain; "
                f"got {self.origin!r} Bounds are {xllim, xrlim} and {yllim, yrlim} respectively"
            )

        return xc, yc

    def _setup_plots(self):
        from matplotlib.mathtext import MathTextParser

        if self._plot_valid:
            return
        if not self._data_valid:
            self._recreate_frb()
        self._colorbar_valid = True
        field_list = list(set(self.data_source._determine_fields(self.fields)))
        for f in field_list:
            axis_index = self.data_source.axis

            xc, yc = self._setup_origin()
            if self.ds._uses_code_length_unit:
                # this should happen only if the dataset was initialized with
                # argument unit_system="code" or if it's set to have no CGS
                # equivalent.  This only needs to happen here in the specific
                # case that we're doing a computationally intense operation
                # like using cartopy, but it prevents crashes in that case.
                (unit_x, unit_y) = ("code_length", "code_length")
            elif self._axes_unit_names is None:
                unit = self.ds.get_smallest_appropriate_unit(
                    self.xlim[1] - self.xlim[0]
                )
                unit_x = unit_y = unit
                coords = self.ds.coordinates
                if hasattr(coords, "image_units"):
                    # check for special cases defined in
                    # non cartesian CoordinateHandler subclasses
                    image_units = coords.image_units[coords.axis_id[axis_index]]
                    if image_units[0] in ("deg", "rad"):
                        unit_x = "code_length"
                    elif image_units[0] == 1:
                        unit_x = "dimensionless"
                    if image_units[1] in ("deg", "rad"):
                        unit_y = "code_length"
                    elif image_units[1] == 1:
                        unit_y = "dimensionless"
            else:
                (unit_x, unit_y) = self._axes_unit_names

            # For some plots we may set aspect by hand, such as for spectral cube data.
            # This will likely be replaced at some point by the coordinate handler
            # setting plot aspect.
            if self.aspect is None:
                self.aspect = float(
                    (self.ds.quan(1.0, unit_y) / self.ds.quan(1.0, unit_x)).in_cgs()
                )
            extentx = (self.xlim - xc)[:2]
            extenty = (self.ylim - yc)[:2]

            # extentx/y arrays inherit units from xlim and ylim attributes
            # and these attributes are always length even for angular and
            # dimensionless axes so we need to stip out units for consistency
            if unit_x == "dimensionless":
                extentx = extentx / extentx.units
            else:
                extentx.convert_to_units(unit_x)
            if unit_y == "dimensionless":
                extenty = extenty / extenty.units
            else:
                extenty.convert_to_units(unit_y)

            extent = [*extentx, *extenty]

            image = self.frb[f]
            font_size = self._font_properties.get_size()

            if f in self.plots.keys():
                pnh = self.plots[f].norm_handler
                cbh = self.plots[f].colorbar_handler
            else:
                pnh, cbh = self._get_default_handlers(
                    field=f, default_display_units=image.units
                )
                if pnh.display_units != image.units:
                    equivalency, equivalency_kwargs = self._equivalencies[f]
                    image.convert_to_units(
                        pnh.display_units, equivalency, **equivalency_kwargs
                    )

            fig = None
            axes = None
            cax = None
            draw_axes = True
            draw_frame = None
            if f in self.plots:
                draw_axes = self.plots[f]._draw_axes
                draw_frame = self.plots[f]._draw_frame
                if self.plots[f].figure is not None:
                    fig = self.plots[f].figure
                    axes = self.plots[f].axes
                    cax = self.plots[f].cax

            # This is for splatting particle positions with a single
            # color instead of a colormap
            if self._splat_color is not None:
                # make image a rgba array, using the splat color
                greyscale_image = self.frb[f]
                ia = np.zeros((greyscale_image.shape[0], greyscale_image.shape[1], 4))
                ia[:, :, 3] = 0.0  # set alpha to 0.0
                locs = greyscale_image > 0.0
                to_rgba = matplotlib.colors.colorConverter.to_rgba
                color_tuple = to_rgba(self._splat_color)
                ia[locs] = color_tuple
                ia = ImageArray(ia)
            else:
                ia = image

            swap_axes = self._has_swapped_axes
            aspect = self.aspect
            if swap_axes:
                extent = _swap_axes_extents(extent)
                ia = ia.transpose()
                aspect = 1.0 / aspect  # aspect ends up passed to imshow(aspect=aspect)

            self.plots[f] = WindowPlotMPL(
                ia,
                extent,
                self.figure_size,
                font_size,
                aspect,
                fig,
                axes,
                cax,
                self._projection,
                self._transform,
                norm_handler=pnh,
                colorbar_handler=cbh,
            )

            axes_unit_labels = self._get_axes_unit_labels(unit_x, unit_y)

            if self.oblique:
                labels = [
                    r"$\rm{Image\ x" + axes_unit_labels[0] + "}$",
                    r"$\rm{Image\ y" + axes_unit_labels[1] + "}$",
                ]
            else:
                coordinates = self.ds.coordinates
                axis_names = coordinates.image_axis_name[axis_index]
                xax = coordinates.x_axis[axis_index]
                yax = coordinates.y_axis[axis_index]

                if hasattr(coordinates, "axis_default_unit_name"):
                    axes_unit_labels = [
                        coordinates.axis_default_unit_name[xax],
                        coordinates.axis_default_unit_name[yax],
                    ]
                labels = [
                    r"$\rm{" + axis_names[0] + axes_unit_labels[0] + r"}$",
                    r"$\rm{" + axis_names[1] + axes_unit_labels[1] + r"}$",
                ]

                if hasattr(coordinates, "axis_field"):
                    if xax in coordinates.axis_field:
                        xmin, xmax = coordinates.axis_field[xax](
                            0, self.xlim, self.ylim
                        )
                    else:
                        xmin, xmax = (float(x) for x in extentx)
                    if yax in coordinates.axis_field:
                        ymin, ymax = coordinates.axis_field[yax](
                            1, self.xlim, self.ylim
                        )
                    else:
                        ymin, ymax = (float(y) for y in extenty)
                    new_extent = (xmin, xmax, ymin, ymax)
                    if swap_axes:
                        new_extent = _swap_axes_extents(new_extent)
                    self.plots[f].image.set_extent(new_extent)
                    self.plots[f].axes.set_aspect("auto")

            x_label, y_label, colorbar_label = self._get_axes_labels(f)

            if x_label is not None:
                labels[0] = x_label
            if y_label is not None:
                labels[1] = y_label

            if swap_axes:
                labels.reverse()

            self.plots[f].axes.set_xlabel(labels[0])
            self.plots[f].axes.set_ylabel(labels[1])

            # Determine the units of the data
            units = Unit(self.frb[f].units, registry=self.ds.unit_registry)
            units = units.latex_representation()

            if colorbar_label is None:
                colorbar_label = image.info["label"]
                if getattr(self, "moment", 1) == 2:
                    colorbar_label = "%s \\rm{Standard Deviation}" % colorbar_label
                if hasattr(self, "projected"):
                    colorbar_label = "$\\rm{Projected }$ %s" % colorbar_label
                if units is None or units == "":
                    pass
                else:
                    colorbar_label += r"$\ \ \left(" + units + r"\right)$"

            parser = MathTextParser("Agg")
            from pyparsing import ParseFatalException

            try:
                parser.parse(colorbar_label)
            except ParseFatalException as err:
                raise YTCannotParseUnitDisplayName(f, colorbar_label, str(err)) from err

            self.plots[f].cb.set_label(colorbar_label)

            # x-y axes minorticks
            if f not in self._minorticks:
                self._minorticks[f] = True
            if self._minorticks[f]:
                self.plots[f].axes.minorticks_on()
            else:
                self.plots[f].axes.minorticks_off()

            if not draw_axes:
                self.plots[f]._toggle_axes(draw_axes, draw_frame)

        self._set_font_properties()
        self.run_callbacks()

        if self._flip_horizontal or self._flip_vertical:
            # some callbacks (e.g., streamlines) fail when applied to a
            # flipped axis, so flip only at the end.
            for f in field_list:
                if self._flip_horizontal:
                    ax = self.plots[f].axes
                    ax.invert_xaxis()

                if self._flip_vertical:
                    ax = self.plots[f].axes
                    ax.invert_yaxis()

        self._plot_valid = True

    def setup_callbacks(self):
        issue_deprecation_warning(
            "The PWViewer.setup_callbacks method is a no-op.", since="4.1.0"
        )

    @invalidate_plot
    def clear_annotations(self, index: Optional[int] = None):
        """
        Clear callbacks from the plot.  If index is not set, clear all
        callbacks.  If index is set, clear that index (ie 0 is the first one
        created, 1 is the 2nd one created, -1 is the last one created, etc.)
        """
        if index is None:
            self._callbacks.clear()
        else:
            self._callbacks.pop(index)
        return self

    def list_annotations(self):
        """
        List the current callbacks for the plot, along with their index.  This
        index can be used with `clear_annotations` to remove a callback from the
        current plot.
        """
        for i, cb in enumerate(self._callbacks):
            print(i, cb)

    def run_callbacks(self):
        for f in self.fields:
            keys = self.frb.keys()
            for callback in self._callbacks:
                # need to pass _swap_axes and adjust all the callbacks
                cbw = CallbackWrapper(
                    self,
                    self.plots[f],
                    self.frb,
                    f,
                    self._font_properties,
                    self._font_color,
                )
                try:
                    callback(cbw)
                except (NotImplementedError, YTDataTypeUnsupported):
                    raise
                except Exception as e:
                    raise YTPlotCallbackError(callback._type_name) from e
            for key in self.frb.keys():
                if key not in keys:
                    del self.frb[key]

    def export_to_mpl_figure(
        self,
        nrows_ncols,
        axes_pad=1.0,
        label_mode="L",
        cbar_location="right",
        cbar_size="5%",
        cbar_mode="each",
        cbar_pad="0%",
    ):
        r"""
        Creates a matplotlib figure object with the specified axes arrangement,
        nrows_ncols, and maps the underlying figures to the matplotlib axes.
        Note that all of these parameters are fed directly to the matplotlib ImageGrid
        class to create the new figure layout.

        Parameters
        ----------

        nrows_ncols : tuple
           the number of rows and columns of the axis grid (e.g., nrows_ncols=(2,2,))
        axes_pad : float
           padding between axes in inches
        label_mode : one of "L", "1", "all"
           arrangement of axes that are labeled
        cbar_location : one of "left", "right", "bottom", "top"
           where to place the colorbar
        cbar_size : string (percentage)
           scaling of the colorbar (e.g., "5%")
        cbar_mode : one of "each", "single", "edge", None
           how to represent the colorbar
        cbar_pad : string (percentage)
           padding between the axis and colorbar (e.g. "5%")

        Returns
        -------

        The return is a matplotlib figure object.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load_sample("IsolatedGalaxy")
        >>> fields = ["density", "velocity_x", "velocity_y", "velocity_magnitude"]
        >>> p = yt.SlicePlot(ds, "z", fields)
        >>> p.set_log("velocity_x", False)
        >>> p.set_log("velocity_y", False)
        >>> fig = p.export_to_mpl_figure((2, 2))
        >>> fig.tight_layout()
        >>> fig.savefig("test.png")

        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import ImageGrid

        fig = plt.figure()
        grid = ImageGrid(
            fig,
            111,
            nrows_ncols=nrows_ncols,
            axes_pad=axes_pad,
            label_mode=label_mode,
            cbar_location=cbar_location,
            cbar_size=cbar_size,
            cbar_mode=cbar_mode,
            cbar_pad=cbar_pad,
        )

        fields = self.fields
        if len(fields) > len(grid):
            raise IndexError("not enough axes for the number of fields")

        for i, f in enumerate(self.fields):
            plot = self.plots[f]
            plot.figure = fig
            plot.axes = grid[i].axes
            plot.cax = grid.cbar_axes[i]

        self._setup_plots()

        return fig


class NormalPlot(abc.ABC):
    """This is the abstraction for SlicePlot and ProjectionPlot, where
    we define the common sanitizing mechanism for user input (normal direction).
    """

    @staticmethod
    def sanitize_normal_vector(ds, normal) -> Union[str, np.ndarray]:
        """Return the name of a cartesian axis whener possible,
        or a 3-element 1D ndarray of float64 in any other valid case.
        Fail with a descriptive error message otherwise.
        """
        axis_names = ds.coordinates.axis_order

        if isinstance(normal, str):
            if normal not in axis_names:
                names_str = ", ".join(f"'{name}'" for name in axis_names)
                raise ValueError(
                    f"'{normal}' is not a valid axis name. Expected one of {names_str}."
                )
            return normal

        if isinstance(normal, int):
            if normal not in (0, 1, 2):
                raise ValueError(
                    f"{normal} is not a valid axis identifier. Expected either 0, 1, or 2."
                )
            return axis_names[normal]

        if not is_sequence(normal):
            raise TypeError(
                f"{normal} is not a valid normal vector identifier. "
                "Expected a string, integer or sequence of 3 floats."
            )

        if len(normal) != 3:
            raise ValueError(
                f"{normal} with length {len(normal)} is not a valid normal vector. "
                "Expected a 3-element sequence."
            )

        try:
            retv = np.array(normal, dtype="float64")
            if retv.shape != (3,):
                raise ValueError(f"{normal} is incorrectly shaped.")
        except ValueError as exc:
            raise TypeError(f"{normal} is not a valid normal vector.") from exc

        nonzero_idx = np.nonzero(retv)[0]
        if len(nonzero_idx) == 0:
            raise ValueError(f"A null vector {normal} isn't a valid normal vector.")
        if len(nonzero_idx) == 1:
            return axis_names[nonzero_idx[0]]

        return retv

    @staticmethod
    def _validate_init_args(*, normal, axis, fields) -> None:
        # TODO: remove this method in yt 4.2

        if axis is not None:
            issue_deprecation_warning(
                "Argument 'axis' is a deprecated alias for 'normal'.",
                since="4.1.0",
                removal="4.2.0",
            )
            if normal is not None:
                raise TypeError("Received incompatible arguments 'axis' and 'normal'")
            normal = axis

        if normal is fields is None:
            raise TypeError(
                "missing 2 required positional arguments: 'normal' and 'fields'"
            )

        if fields is None:
            raise TypeError("missing required positional argument: 'fields'")

        if normal is None:
            raise TypeError("missing required positional argument: 'normal'")

        return normal


class SlicePlot(NormalPlot):
    r"""
    A dispatch class for :class:`yt.visualization.plot_window.AxisAlignedSlicePlot`
    and :class:`yt.visualization.plot_window.OffAxisSlicePlot` objects.  This
    essentially allows for a single entry point to both types of slice plots,
    the distinction being determined by the specified normal vector to the
    projection.

    The returned plot object can be updated using one of the many helper
    functions defined in PlotWindow.

    Parameters
    ----------

    ds : :class:`yt.data_objects.static_output.Dataset`
        This is the dataset object corresponding to the
        simulation output to be plotted.
    normal : int, str, or 3-element sequence of floats
        This specifies the normal vector to the slice.
        Valid int values are 0, 1 and 2. Coresponding str values depend on the
        geometry of the dataset and are generally given by `ds.coordinates.axis_order`.
        E.g. in cartesian they are 'x', 'y' and 'z'.
        An arbitrary normal vector may be specified as a 3-element sequence of floats.

        This returns a :class:`OffAxisSlicePlot` object or a
        :class:`AxisAlignedSlicePlot` object, depending on wether the requested
        normal directions corresponds to a natural axis of the dataset's geometry.

    fields : a (or a list of) 2-tuple of strings (ftype, fname)
         The name of the field(s) to be plotted.

    The following are nominally keyword arguments passed onto the respective
    slice plot objects generated by this function.

    Keyword Arguments
    -----------------

    center : A sequence floats, a string, or a tuple.
         The coordinate of the center of the image. If set to 'c', 'center' or
         left blank, the plot is centered on the middle of the domain. If set to
         'max' or 'm', the center will be located at the maximum of the
         ('gas', 'density') field. Centering on the max or min of a specific
         field is supported by providing a tuple such as ("min","temperature") or
         ("max","dark_matter_density"). Units can be specified by passing in *center*
         as a tuple containing a coordinate and string unit name or by passing
         in a YTArray. If a list or unitless array is supplied, code units are
         assumed.
    width : tuple or a float.
         Width can have four different formats to support windows with variable
         x and y widths.  They are:

         ==================================     =======================
         format                                 example
         ==================================     =======================
         (float, string)                        (10,'kpc')
         ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
         float                                  0.2
         (float, float)                         (0.2, 0.3)
         ==================================     =======================

         For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
         wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
         window that is 10 kiloparsecs wide along the x axis and 15
         kiloparsecs wide along the y axis.  In the other two examples, code
         units are assumed, for example (0.2, 0.3) requests a plot that has an
         x width of 0.2 and a y width of 0.3 in code units.  If units are
         provided the resulting plot axis labels will use the supplied units.
    axes_unit : string
         The name of the unit for the tick labels on the x and y axes.
         Defaults to None, which automatically picks an appropriate unit.
         If axes_unit is '1', 'u', or 'unitary', it will not display the
         units, and only show the axes name.
    origin : string or length 1, 2, or 3 sequence.
         The location of the origin of the plot coordinate system for
         `AxisAlignedSlicePlot` object; for `OffAxisSlicePlot` objects this
         parameter is discarded. This is typically represented by a '-'
         separated string or a tuple of strings. In the first index the
         y-location is given by 'lower', 'upper', or 'center'. The second index
         is the x-location, given as 'left', 'right', or 'center'. Finally, the
         whether the origin is applied in 'domain' space, plot 'window' space or
         'native' simulation coordinate system is given. For example, both
         'upper-right-domain' and ['upper', 'right', 'domain'] place the
         origin in the upper right hand corner of domain space. If x or y
         are not given, a value is inferred. For instance, 'left-domain'
         corresponds to the lower-left hand corner of the simulation domain,
         'center-domain' corresponds to the center of the simulation domain,
         or 'center-window' for the center of the plot window. In the event
         that none of these options place the origin in a desired location,
         a sequence of tuples and a string specifying the
         coordinate space can be given. If plain numeric types are input,
         units of `code_length` are assumed. Further examples:

         =============================================== ===============================
         format                                          example
         =============================================== ===============================
         '{space}'                                       'domain'
         '{xloc}-{space}'                                'left-window'
         '{yloc}-{space}'                                'upper-domain'
         '{yloc}-{xloc}-{space}'                         'lower-right-window'
         ('{space}',)                                    ('window',)
         ('{xloc}', '{space}')                           ('right', 'domain')
         ('{yloc}', '{space}')                           ('lower', 'window')
         ('{yloc}', '{xloc}', '{space}')                 ('lower', 'right', 'window')
         ((yloc, '{unit}'), (xloc, '{unit}'), '{space}') ((0, 'm'), (.4, 'm'), 'window')
         (xloc, yloc, '{space}')                         (0.23, 0.5, 'domain')
         =============================================== ===============================
    north_vector : a sequence of floats
        A vector defining the 'up' direction in the `OffAxisSlicePlot`; not
        used in `AxisAlignedSlicePlot`.  This option sets the orientation of the
        slicing plane.  If not set, an arbitrary grid-aligned north-vector is
        chosen.
    fontsize : integer
         The size of the fonts for the axis, colorbar, and tick labels.
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source : YTSelectionContainer Object
         Object to be used for data selection.  Defaults to a region covering
         the entire simulation.
    swap_axes : bool


    Raises
    ------

    ValueError or TypeError
        If `normal` cannot be interpreted as a valid normal direction.

    Examples
    --------

    >>> from yt import load
    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> slc = SlicePlot(ds, "x", ("gas", "density"), center=[0.2, 0.3, 0.4])

    >>> slc = SlicePlot(
    ...     ds, [0.4, 0.2, -0.1], ("gas", "pressure"), north_vector=[0.2, -0.3, 0.1]
    ... )

    """

    # ignoring type check here, because mypy doesn't allow __new__ methods to
    # return instances of subclasses. The design we use here is however based
    # on the pathlib.Path class from the standard library
    # https://github.com/python/mypy/issues/1020
    def __new__(  # type: ignore
        cls, ds, normal=None, fields=None, *args, axis=None, **kwargs
    ) -> Union["AxisAlignedSlicePlot", "OffAxisSlicePlot"]:
        # TODO: in yt 4.2, remove default values for normal and fields, drop axis kwarg
        normal = cls._validate_init_args(normal=normal, axis=axis, fields=fields)

        if cls is SlicePlot:
            normal = cls.sanitize_normal_vector(ds, normal)
            if isinstance(normal, str):
                cls = AxisAlignedSlicePlot
            else:
                cls = OffAxisSlicePlot
        self = object.__new__(cls)
        return self  # type: ignore [return-value]


class ProjectionPlot(NormalPlot):
    r"""
    A dispatch class for :class:`yt.visualization.plot_window.AxisAlignedProjectionPlot`
    and :class:`yt.visualization.plot_window.OffAxisProjectionPlot` objects.  This
    essentially allows for a single entry point to both types of projection plots,
    the distinction being determined by the specified normal vector to the
    slice.

    The returned plot object can be updated using one of the many helper
    functions defined in PlotWindow.

    Parameters
    ----------

    ds : :class:`yt.data_objects.static_output.Dataset`
        This is the dataset object corresponding to the
        simulation output to be plotted.
    normal : int, str, or 3-element sequence of floats
        This specifies the normal vector to the slice.
        Valid int values are 0, 1 and 2. Coresponding str values depend on the
        geometry of the dataset and are generally given by `ds.coordinates.axis_order`.
        E.g. in cartesian they are 'x', 'y' and 'z'.
        An arbitrary normal vector may be specified as a 3-element sequence of floats.

        This function will return a :class:`OffAxisProjectionPlot` object or a
        :class:`AxisAlignedProjectionPlot` object, depending on wether the requested
        normal directions corresponds to a natural axis of the dataset's geometry.

    fields : a (or a list of) 2-tuple of strings (ftype, fname)
         The name of the field(s) to be plotted.


    Any additional positional and keyword arguments are passed down to the appropriate
    return class. See :class:`yt.visualization.plot_window.AxisAlignedProjectionPlot`
    and :class:`yt.visualization.plot_window.OffAxisProjectionPlot`.

    Raises
    ------

    ValueError or TypeError
        If `normal` cannot be interpreted as a valid normal direction.

    """

    # ignoring type check here, because mypy doesn't allow __new__ methods to
    # return instances of subclasses. The design we use here is however based
    # on the pathlib.Path class from the standard library
    # https://github.com/python/mypy/issues/1020
    def __new__(  # type: ignore
        cls, ds, normal=None, fields=None, *args, axis=None, **kwargs
    ) -> Union["AxisAlignedProjectionPlot", "OffAxisProjectionPlot"]:
        # TODO: in yt 4.2, remove default values for normal and fields, drop axis kwarg
        normal = cls._validate_init_args(normal=normal, axis=axis, fields=fields)

        if cls is ProjectionPlot:
            normal = cls.sanitize_normal_vector(ds, normal)
            if isinstance(normal, str):
                cls = AxisAlignedProjectionPlot
            else:
                cls = OffAxisProjectionPlot
        self = object.__new__(cls)
        return self  # type: ignore [return-value]


class AxisAlignedSlicePlot(SlicePlot, PWViewerMPL):
    r"""Creates a slice plot from a dataset

    Given a ds object, an axis to slice along, and a field name
    string, this will return a PWViewerMPL object containing
    the plot.

    The plot can be updated using one of the many helper functions
    defined in PlotWindow.

    Parameters
    ----------
    ds : `Dataset`
         This is the dataset object corresponding to the
         simulation output to be plotted.
    normal : int or one of 'x', 'y', 'z'
         An int corresponding to the axis to slice along (0=x, 1=y, 2=z)
         or the axis name itself
    fields : string
         The name of the field(s) to be plotted.
    center : A sequence of floats, a string, or a tuple.
         The coordinate of the center of the image. If set to 'c', 'center' or
         left blank, the plot is centered on the middle of the domain. If set to
         'max' or 'm', the center will be located at the maximum of the
         ('gas', 'density') field. Centering on the max or min of a specific
         field is supported by providing a tuple such as ("min","temperature") or
         ("max","dark_matter_density"). Units can be specified by passing in *center*
         as a tuple containing a coordinate and string unit name or by passing
         in a YTArray. If a list or unitless array is supplied, code units are
         assumed.
    width : tuple or a float.
         Width can have four different formats to support windows with variable
         x and y widths.  They are:

         ==================================     =======================
         format                                 example
         ==================================     =======================
         (float, string)                        (10,'kpc')
         ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
         float                                  0.2
         (float, float)                         (0.2, 0.3)
         ==================================     =======================

         For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
         wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
         window that is 10 kiloparsecs wide along the x axis and 15
         kiloparsecs wide along the y axis.  In the other two examples, code
         units are assumed, for example (0.2, 0.3) requests a plot that has an
         x width of 0.2 and a y width of 0.3 in code units.  If units are
         provided the resulting plot axis labels will use the supplied units.
    origin : string or length 1, 2, or 3 sequence.
         The location of the origin of the plot coordinate system. This
         is typically represented by a '-' separated string or a tuple of
         strings. In the first index the y-location is given by 'lower',
         'upper', or 'center'. The second index is the x-location, given as
         'left', 'right', or 'center'. Finally, whether the origin is
         applied in 'domain' space, plot 'window' space or 'native'
         simulation coordinate system is given. For example, both
         'upper-right-domain' and ['upper', 'right', 'domain'] place the
         origin in the upper right hand corner of domain space. If x or y
         are not given, a value is inferred. For instance, 'left-domain'
         corresponds to the lower-left hand corner of the simulation domain,
         'center-domain' corresponds to the center of the simulation domain,
         or 'center-window' for the center of the plot window. In the event
         that none of these options place the origin in a desired location,
         a sequence of tuples and a string specifying the
         coordinate space can be given. If plain numeric types are input,
         units of `code_length` are assumed. Further examples:

         =============================================== ===============================
         format                                          example
         =============================================== ===============================
         '{space}'                                       'domain'
         '{xloc}-{space}'                                'left-window'
         '{yloc}-{space}'                                'upper-domain'
         '{yloc}-{xloc}-{space}'                         'lower-right-window'
         ('{space}',)                                    ('window',)
         ('{xloc}', '{space}')                           ('right', 'domain')
         ('{yloc}', '{space}')                           ('lower', 'window')
         ('{yloc}', '{xloc}', '{space}')                 ('lower', 'right', 'window')
         ((yloc, '{unit}'), (xloc, '{unit}'), '{space}') ((0, 'm'), (.4, 'm'), 'window')
         (xloc, yloc, '{space}')                         (0.23, 0.5, 'domain')
         =============================================== ===============================
    axes_unit : string
         The name of the unit for the tick labels on the x and y axes.
         Defaults to None, which automatically picks an appropriate unit.
         If axes_unit is '1', 'u', or 'unitary', it will not display the
         units, and only show the axes name.
    right_handed : boolean
        Depreceated, please use flip_horizontal callback.
        Whether the implicit east vector for the image generated is set to make a right
        handed coordinate system with a north vector and the normal vector, the
        direction of the 'window' into the data.
    fontsize : integer
         The size of the fonts for the axis, colorbar, and tick labels.
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source: YTSelectionContainer object
         Object to be used for data selection. Defaults to ds.all_data(), a
         region covering the full domain
    buff_size: length 2 sequence
         Size of the buffer to use for the image, i.e. the number of resolution elements
         used.  Effectively sets a resolution limit to the image if buff_size is
         smaller than the finest gridding.

    Examples
    --------

    This will save an image in the file 'sliceplot_Density.png'

    >>> from yt import load
    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> p = SlicePlot(ds, 2, "density", "c", (20, "kpc"))
    >>> p.save("sliceplot")

    """
    _plot_type = "Slice"
    _frb_generator = FixedResolutionBuffer

    def __init__(
        self,
        ds,
        normal=None,
        fields=None,
        center="c",
        width=None,
        axes_unit=None,
        origin="center-window",
        right_handed=True,
        fontsize=18,
        field_parameters=None,
        window_size=8.0,
        aspect=None,
        data_source=None,
        buff_size=(800, 800),
        *,
        north_vector=None,
        axis=None,
    ):
        # TODO: in yt 4.2, remove default values for normal and fields, drop axis kwarg
        if north_vector is not None:
            # this kwarg exists only for symmetry reasons with OffAxisSlicePlot
            mylog.warning(
                "Ignoring 'north_vector' keyword as it is ill-defined for "
                "an AxisAlignedSlicePlot object."
            )
            del north_vector

        normal = self._validate_init_args(
            normal=normal,
            axis=axis,
            fields=fields,
        )
        normal = self.sanitize_normal_vector(ds, normal)
        # this will handle time series data and controllers
        axis = fix_axis(normal, ds)
        (bounds, center, display_center) = get_window_parameters(
            axis, center, width, ds
        )
        if field_parameters is None:
            field_parameters = {}

        if isinstance(ds, YTSpatialPlotDataset):
            slc = ds.all_data()
            slc.axis = axis
            if slc.axis != ds.parameters["axis"]:
                raise RuntimeError(f"Original slice axis is {ds.parameters['axis']}.")
        else:
            slc = ds.slice(
                axis,
                center[axis],
                field_parameters=field_parameters,
                center=center,
                data_source=data_source,
            )
            slc.get_data(fields)
        validate_mesh_fields(slc, fields)
        PWViewerMPL.__init__(
            self,
            slc,
            bounds,
            origin=origin,
            fontsize=fontsize,
            fields=fields,
            window_size=window_size,
            aspect=aspect,
            right_handed=right_handed,
            buff_size=buff_size,
            geometry=ds.geometry,
        )
        if axes_unit is None:
            axes_unit = get_axes_unit(width, ds)
        self.set_axes_unit(axes_unit)


class AxisAlignedProjectionPlot(ProjectionPlot, PWViewerMPL):
    r"""Creates a projection plot from a dataset

    Given a ds object, an axis to project along, and a field name
    string, this will return a PWViewerMPL object containing
    the plot.

    The plot can be updated using one of the many helper functions
    defined in PlotWindow.

    Parameters
    ----------
    ds : `Dataset`
        This is the dataset object corresponding to the
        simulation output to be plotted.
    normal : int or one of 'x', 'y', 'z'
        An int corresponding to the axis to slice along (0=x, 1=y, 2=z)
        or the axis name itself
    fields : string
        The name of the field(s) to be plotted.
    center : A sequence of floats, a string, or a tuple.
        The coordinate of the center of the image. If set to 'c', 'center' or
        left blank, the plot is centered on the middle of the domain. If set to
        'max' or 'm', the center will be located at the maximum of the
        ('gas', 'density') field. Centering on the max or min of a specific
        field is supported by providing a tuple such as ("min","temperature") or
        ("max","dark_matter_density"). Units can be specified by passing in *center*
        as a tuple containing a coordinate and string unit name or by passing
        in a YTArray. If a list or unitless array is supplied, code units are
        assumed.
    width : tuple or a float.
        Width can have four different formats to support windows with variable
        x and y widths.  They are:

        ==================================     =======================
        format                                 example
        ==================================     =======================
        (float, string)                        (10,'kpc')
        ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
        float                                  0.2
        (float, float)                         (0.2, 0.3)
        ==================================     =======================

        For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
        wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
        window that is 10 kiloparsecs wide along the x axis and 15
        kiloparsecs wide along the y axis.  In the other two examples, code
        units are assumed, for example (0.2, 0.3) requests a plot that has an
        x width of 0.2 and a y width of 0.3 in code units.  If units are
        provided the resulting plot axis labels will use the supplied units.
    axes_unit : string
        The name of the unit for the tick labels on the x and y axes.
        Defaults to None, which automatically picks an appropriate unit.
        If axes_unit is '1', 'u', or 'unitary', it will not display the
        units, and only show the axes name.
    origin : string or length 1, 2, or 3 sequence.
        The location of the origin of the plot coordinate system. This
        is typically represented by a '-' separated string or a tuple of
        strings. In the first index the y-location is given by 'lower',
        'upper', or 'center'. The second index is the x-location, given as
        'left', 'right', or 'center'. Finally, whether the origin is
        applied in 'domain' space, plot 'window' space or 'native'
        simulation coordinate system is given. For example, both
        'upper-right-domain' and ['upper', 'right', 'domain'] place the
        origin in the upper right hand corner of domain space. If x or y
        are not given, a value is inferred. For instance, 'left-domain'
        corresponds to the lower-left hand corner of the simulation domain,
        'center-domain' corresponds to the center of the simulation domain,
        or 'center-window' for the center of the plot window. In the event
        that none of these options place the origin in a desired location,
        a sequence of tuples and a string specifying the
        coordinate space can be given. If plain numeric types are input,
        units of `code_length` are assumed. Further examples:

        =============================================== ===============================
        format                                          example
        =============================================== ===============================
        '{space}'                                       'domain'
        '{xloc}-{space}'                                'left-window'
        '{yloc}-{space}'                                'upper-domain'
        '{yloc}-{xloc}-{space}'                         'lower-right-window'
        ('{space}',)                                    ('window',)
        ('{xloc}', '{space}')                           ('right', 'domain')
        ('{yloc}', '{space}')                           ('lower', 'window')
        ('{yloc}', '{xloc}', '{space}')                 ('lower', 'right', 'window')
        ((yloc, '{unit}'), (xloc, '{unit}'), '{space}') ((0, 'm'), (.4, 'm'), 'window')
        (xloc, yloc, '{space}')                            (0.23, 0.5, 'domain')
        =============================================== ===============================

    right_handed : boolean
        Depreceated, please use flip_horizontal callback.
        Whether the implicit east vector for the image generated is set to make a right
        handed coordinate system with a north vector and the normal vector, the
        direction of the 'window' into the data.
    data_source : YTSelectionContainer Object
        Object to be used for data selection.  Defaults to a region covering
        the entire simulation.
    weight_field : string
        The name of the weighting field.  Set to None for no weight.
    max_level: int
        The maximum level to project to.
    fontsize : integer
        The size of the fonts for the axis, colorbar, and tick labels.
    method : string
        The method of projection.  Valid methods are:

        "integrate" with no weight_field specified : integrate the requested
        field along the line of sight.

        "integrate" with a weight_field specified : weight the requested
        field by the weighting field and integrate along the line of sight.

        "max" : pick out the maximum value of the field in the line of sight.
        "min" : pick out the minimum value of the field in the line of sight.

        "sum" : This method is the same as integrate, except that it does not
        multiply by a path length when performing the integration, and is
        just a straight summation of the field along the given axis. WARNING:
        This should only be used for uniform resolution grid datasets, as other
        datasets may result in unphysical images.
    window_size : float
        The size of the window in inches. Set to 8 by default.
    aspect : float
        The aspect ratio of the plot.  Set to None for 1.
    field_parameters : dictionary
        A dictionary of field parameters than can be accessed by derived
        fields.
    data_source: YTSelectionContainer object
        Object to be used for data selection. Defaults to ds.all_data(), a
        region covering the full domain
    buff_size: length 2 sequence
        Size of the buffer to use for the image, i.e. the number of resolution elements
        used. Effectively sets a resolution limit to the image if buff_size is
        smaller than the finest gridding.
    moment : integer, optional
        for a weighted projection, moment = 1 (the default) corresponds to a
        weighted average. moment = 2 corresponds to a weighted standard
        deviation.

    Examples
    --------

    Create a projection plot with a width of 20 kiloparsecs centered on the
    center of the simulation box:

    >>> from yt import load
    >>> ds = load("IsolateGalaxygalaxy0030/galaxy0030")
    >>> p = AxisAlignedProjectionPlot(ds, "z", ("gas", "density"), width=(20, "kpc"))

    """
    _plot_type = "Projection"
    _frb_generator = FixedResolutionBuffer

    def __init__(
        self,
        ds,
        normal=None,
        fields=None,
        center="c",
        width=None,
        axes_unit=None,
        weight_field=None,
        max_level=None,
        origin="center-window",
        right_handed=True,
        fontsize=18,
        field_parameters=None,
        data_source=None,
        method="integrate",
        window_size=8.0,
        buff_size=(800, 800),
        aspect=None,
        *,
        moment=1,
        axis=None,
    ):
        if method == "mip":
            issue_deprecation_warning(
                "'mip' method is a deprecated alias for 'max'. "
                "Please use method='max' directly.",
                since="4.1.0",
            )
            method = "max"
        # TODO: in yt 4.2, remove default values for normal and fields, drop axis kwarg
        normal = self._validate_init_args(normal=normal, fields=fields, axis=axis)
        normal = self.sanitize_normal_vector(ds, normal)

        axis = fix_axis(normal, ds)
        # If a non-weighted integral projection, assure field-label reflects that
        if weight_field is None and method == "integrate":
            self.projected = True
        (bounds, center, display_center) = get_window_parameters(
            axis, center, width, ds
        )
        if field_parameters is None:
            field_parameters = {}

        # We don't use the plot's data source for validation like in the other
        # plotting classes to avoid an exception
        test_data_source = ds.all_data()
        validate_mesh_fields(test_data_source, fields)

        if isinstance(ds, YTSpatialPlotDataset):
            proj = ds.all_data()
            proj.axis = axis
            if proj.axis != ds.parameters["axis"]:
                raise RuntimeError(
                    f"Original projection axis is {ds.parameters['axis']}."
                )
            if weight_field is not None:
                proj.weight_field = proj._determine_fields(weight_field)[0]
            else:
                proj.weight_field = weight_field
            proj.center = center
        else:
            proj = ds.proj(
                fields,
                axis,
                weight_field=weight_field,
                center=center,
                data_source=data_source,
                field_parameters=field_parameters,
                method=method,
                max_level=max_level,
                moment=moment,
            )
        self.moment = moment
        PWViewerMPL.__init__(
            self,
            proj,
            bounds,
            fields=fields,
            origin=origin,
            right_handed=right_handed,
            fontsize=fontsize,
            window_size=window_size,
            aspect=aspect,
            buff_size=buff_size,
            geometry=ds.geometry,
        )
        if axes_unit is None:
            axes_unit = get_axes_unit(width, ds)
        self.set_axes_unit(axes_unit)


class OffAxisSlicePlot(SlicePlot, PWViewerMPL):
    r"""Creates an off axis slice plot from a dataset

    Given a ds object, a normal vector defining a slicing plane, and
    a field name string, this will return a PWViewerMPL object
    containing the plot.

    The plot can be updated using one of the many helper functions
    defined in PlotWindow.

    Parameters
    ----------
    ds : :class:`yt.data_objects.static_output.Dataset`
         This is the dataset object corresponding to the
         simulation output to be plotted.
    normal : a sequence of floats
         The vector normal to the slicing plane.
    fields : string
         The name of the field(s) to be plotted.
    center : A sequence of floats, a string, or a tuple.
         The coordinate of the center of the image. If set to 'c', 'center' or
         left blank, the plot is centered on the middle of the domain. If set to
         'max' or 'm', the center will be located at the maximum of the
         ('gas', 'density') field. Centering on the max or min of a specific
         field is supported by providing a tuple such as ("min","temperature") or
         ("max","dark_matter_density"). Units can be specified by passing in *center*
         as a tuple containing a coordinate and string unit name or by passing
         in a YTArray. If a list or unitless array is supplied, code units are
         assumed.
    width : tuple or a float.
         Width can have four different formats to support windows with variable
         x and y widths.  They are:

         ==================================     =======================
         format                                 example
         ==================================     =======================
         (float, string)                        (10,'kpc')
         ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
         float                                  0.2
         (float, float)                         (0.2, 0.3)
         ==================================     =======================

         For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
         wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
         window that is 10 kiloparsecs wide along the x axis and 15
         kiloparsecs wide along the y axis.  In the other two examples, code
         units are assumed, for example (0.2, 0.3) requests a plot that has an
         x width of 0.2 and a y width of 0.3 in code units.  If units are
         provided the resulting plot axis labels will use the supplied units.
    axes_unit : string
         The name of the unit for the tick labels on the x and y axes.
         Defaults to None, which automatically picks an appropriate unit.
         If axes_unit is '1', 'u', or 'unitary', it will not display the
         units, and only show the axes name.
    north_vector : a sequence of floats
         A vector defining the 'up' direction in the plot.  This
         option sets the orientation of the slicing plane.  If not
         set, an arbitrary grid-aligned north-vector is chosen.
    right_handed : boolean
        Depreceated, please use flip_horizontal callback.
        Whether the implicit east vector for the image generated is set to make a right
        handed coordinate system with a north vector and the normal vector, the
        direction of the 'window' into the data.
    fontsize : integer
         The size of the fonts for the axis, colorbar, and tick labels.
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source : YTSelectionContainer Object
         Object to be used for data selection.  Defaults ds.all_data(), a
         region covering the full domain.
    buff_size: length 2 sequence
         Size of the buffer to use for the image, i.e. the number of resolution elements
         used.  Effectively sets a resolution limit to the image if buff_size is
         smaller than the finest gridding.
    """

    _plot_type = "OffAxisSlice"
    _frb_generator = FixedResolutionBuffer
    _supported_geometries = ("cartesian", "spectral_cube")

    def __init__(
        self,
        ds,
        normal,
        fields,
        center="c",
        width=None,
        axes_unit=None,
        north_vector=None,
        right_handed=True,
        fontsize=18,
        field_parameters=None,
        data_source=None,
        buff_size=(800, 800),
        *,
        origin=None,
    ):
        if origin is not None:
            # this kwarg exists only for symmetry reasons with AxisAlignedSlicePlot
            # in OffAxisSlicePlot, the origin is hardcoded
            mylog.warning(
                "Ignoring 'origin' keyword as it is ill-defined for "
                "an OffAxisSlicePlot object."
            )
            del origin

        if ds.geometry not in self._supported_geometries:
            raise NotImplementedError(
                f"off-axis slices are not supported for {ds.geometry!r} geometry\n"
                f"currently supported geometries: {self._supported_geometries!r}"
            )

        (bounds, center_rot) = get_oblique_window_parameters(normal, center, width, ds)
        if field_parameters is None:
            field_parameters = {}

        if isinstance(ds, YTSpatialPlotDataset):
            cutting = ds.all_data()
            cutting.axis = 4
            cutting._inv_mat = ds.parameters["_inv_mat"]
        else:
            cutting = ds.cutting(
                normal,
                center,
                north_vector=north_vector,
                field_parameters=field_parameters,
                data_source=data_source,
            )
            cutting.get_data(fields)
        validate_mesh_fields(cutting, fields)
        # Hard-coding the origin keyword since the other two options
        # aren't well-defined for off-axis data objects
        PWViewerMPL.__init__(
            self,
            cutting,
            bounds,
            fields=fields,
            origin="center-window",
            periodic=False,
            right_handed=right_handed,
            oblique=True,
            fontsize=fontsize,
            buff_size=buff_size,
        )
        if axes_unit is None:
            axes_unit = get_axes_unit(width, ds)
        self.set_axes_unit(axes_unit)


class OffAxisProjectionDummyDataSource:
    _type_name = "proj"
    _key_fields: List[str] = []

    def __init__(
        self,
        center,
        ds,
        normal_vector,
        width,
        fields,
        interpolated,
        weight=None,
        volume=None,
        no_ghost=False,
        le=None,
        re=None,
        north_vector=None,
        method="integrate",
        data_source=None,
        *,
        moment=1,
    ):
        validate_moment(moment, weight)
        self.center = center
        self.ds = ds
        self.axis = 4  # always true for oblique data objects
        self.normal_vector = normal_vector
        self.width = width
        if data_source is None:
            self.dd = ds.all_data()
        else:
            self.dd = data_source
        fields = self.dd._determine_fields(fields)
        self.fields = fields
        self.interpolated = interpolated
        if weight is not None:
            weight = self.dd._determine_fields(weight)[0]
        self.weight_field = weight
        self.volume = volume
        self.no_ghost = no_ghost
        self.le = le
        self.re = re
        self.north_vector = north_vector
        self.method = method
        self.orienter = Orientation(normal_vector, north_vector=north_vector)
        self.moment = moment

    def _determine_fields(self, *args):
        return self.dd._determine_fields(*args)


class OffAxisProjectionPlot(ProjectionPlot, PWViewerMPL):
    r"""Creates an off axis projection plot from a dataset

    Given a ds object, a normal vector to project along, and
    a field name string, this will return a PWViewerMPL object
    containing the plot.

    The plot can be updated using one of the many helper functions
    defined in PlotWindow.

    Parameters
    ----------
    ds : :class:`yt.data_objects.static_output.Dataset`
        This is the dataset object corresponding to the
        simulation output to be plotted.
    normal : a sequence of floats
        The vector normal to the slicing plane.
    fields : string
        The name of the field(s) to be plotted.
    center : A sequence of floats, a string, or a tuple.
        The coordinate of the center of the image. If set to 'c', 'center' or
        left blank, the plot is centered on the middle of the domain. If set to
        'max' or 'm', the center will be located at the maximum of the
        ('gas', 'density') field. Centering on the max or min of a specific
        field is supported by providing a tuple such as ("min","temperature") or
        ("max","dark_matter_density"). Units can be specified by passing in *center*
        as a tuple containing a coordinate and string unit name or by passing
        in a YTArray. If a list or unitless array is supplied, code units are
        assumed.
    width : tuple or a float.
        Width can have four different formats to support windows with variable
        x and y widths. They are:

        ==================================     =======================
        format                                 example
        ==================================     =======================
        (float, string)                        (10,'kpc')
        ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
        float                                  0.2
        (float, float)                         (0.2, 0.3)
        ==================================     =======================

        For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
        wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
        window that is 10 kiloparsecs wide along the x axis and 15
        kiloparsecs wide along the y axis.  In the other two examples, code
        units are assumed, for example (0.2, 0.3) requests a plot that has an
        x width of 0.2 and a y width of 0.3 in code units.  If units are
        provided the resulting plot axis labels will use the supplied units.
    depth : A tuple or a float
        A tuple containing the depth to project through and the string
        key of the unit: (width, 'unit'). If set to a float, code units
        are assumed
    weight_field : string
        The name of the weighting field.  Set to None for no weight.
    max_level: int
        The maximum level to project to.
    axes_unit : string
        The name of the unit for the tick labels on the x and y axes.
        Defaults to None, which automatically picks an appropriate unit.
        If axes_unit is '1', 'u', or 'unitary', it will not display the
        units, and only show the axes name.
    north_vector : a sequence of floats
        A vector defining the 'up' direction in the plot. This
        option sets the orientation of the slicing plane. If not
        set, an arbitrary grid-aligned north-vector is chosen.
    right_handed : boolean
        Depreceated, please use flip_horizontal callback.
        Whether the implicit east vector for the image generated is set to make a right
        handed coordinate system with a north vector and the normal vector, the
        direction of the 'window' into the data.
    fontsize : integer
        The size of the fonts for the axis, colorbar, and tick labels.
    method : string
        The method of projection. Valid methods are:

        "integrate" with no weight_field specified : integrate the requested
        field along the line of sight.

        "integrate" with a weight_field specified : weight the requested
        field by the weighting field and integrate along the line of sight.

        "sum" : This method is the same as integrate, except that it does not
        multiply by a path length when performing the integration, and is
        just a straight summation of the field along the given axis. WARNING:
        This should only be used for uniform resolution grid datasets, as other
        datasets may result in unphysical images.
    moment : integer, optional
        for a weighted projection, moment = 1 (the default) corresponds to a
        weighted average. moment = 2 corresponds to a weighted standard
        deviation.
    data_source: YTSelectionContainer object
        Object to be used for data selection. Defaults to ds.all_data(), a
        region covering the full domain
    buff_size: length 2 sequence
        Size of the buffer to use for the image, i.e. the number of resolution elements
        used. Effectively sets a resolution limit to the image if buff_size is
        smaller than the finest gridding.
    """
    _plot_type = "OffAxisProjection"
    _frb_generator = OffAxisProjectionFixedResolutionBuffer
    _supported_geometries = ("cartesian", "spectral_cube")

    def __init__(
        self,
        ds,
        normal,
        fields,
        center="c",
        width=None,
        depth=(1, "1"),
        axes_unit=None,
        weight_field=None,
        max_level=None,
        north_vector=None,
        right_handed=True,
        volume=None,
        no_ghost=False,
        le=None,
        re=None,
        interpolated=False,
        fontsize=18,
        method="integrate",
        moment=1,
        data_source=None,
        buff_size=(800, 800),
    ):
        if ds.geometry not in self._supported_geometries:
            raise NotImplementedError(
                f"off-axis slices are not supported for {ds.geometry!r} geometry\n"
                f"currently supported geometries: {self._supported_geometries!r}"
            )

        (bounds, center_rot) = get_oblique_window_parameters(
            normal, center, width, ds, depth=depth
        )
        fields = list(iter_fields(fields))[:]
        oap_width = ds.arr(
            (bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
        )
        OffAxisProj = OffAxisProjectionDummyDataSource(
            center_rot,
            ds,
            normal,
            oap_width,
            fields,
            interpolated,
            weight=weight_field,
            volume=volume,
            no_ghost=no_ghost,
            le=le,
            re=re,
            north_vector=north_vector,
            method=method,
            data_source=data_source,
            moment=moment,
        )

        validate_mesh_fields(OffAxisProj, fields)

        if max_level is not None:
            OffAxisProj.dd.max_level = max_level

        # If a non-weighted, integral projection, assure field label
        # reflects that
        if weight_field is None and OffAxisProj.method == "integrate":
            self.projected = True

        self.moment = moment

        # Hard-coding the origin keyword since the other two options
        # aren't well-defined for off-axis data objects
        PWViewerMPL.__init__(
            self,
            OffAxisProj,
            bounds,
            fields=fields,
            origin="center-window",
            periodic=False,
            oblique=True,
            right_handed=right_handed,
            fontsize=fontsize,
            buff_size=buff_size,
        )
        if axes_unit is None:
            axes_unit = get_axes_unit(width, ds)
        self.set_axes_unit(axes_unit)


class WindowPlotMPL(ImagePlotMPL):
    """A container for a single PlotWindow matplotlib figure and axes"""

    def __init__(
        self,
        data,
        extent,
        figure_size,
        fontsize,
        aspect,
        figure,
        axes,
        cax,
        mpl_proj,
        mpl_transform,
        *,
        norm_handler: NormHandler,
        colorbar_handler: ColorbarHandler,
    ):
        self._projection = mpl_proj
        self._transform = mpl_transform

        self._setup_layout_constraints(figure_size, fontsize)
        self._draw_frame = True
        self._aspect = ((extent[1] - extent[0]) / (extent[3] - extent[2])).in_cgs()
        self._unit_aspect = aspect

        # Compute layout
        self._figure_size = figure_size
        self._draw_axes = True
        fontscale = float(fontsize) / self.__class__._default_font_size
        if fontscale < 1.0:
            fontscale = np.sqrt(fontscale)

        if is_sequence(figure_size):
            self._cb_size = 0.0375 * figure_size[0]
        else:
            self._cb_size = 0.0375 * figure_size
        self._ax_text_size = [1.2 * fontscale, 0.9 * fontscale]
        self._top_buff_size = 0.30 * fontscale

        super().__init__(
            figure=figure,
            axes=axes,
            cax=cax,
            norm_handler=norm_handler,
            colorbar_handler=colorbar_handler,
        )

        self._init_image(data, extent, aspect)

    def _create_axes(self, axrect):
        self.axes = self.figure.add_axes(axrect, projection=self._projection)


def plot_2d(
    ds,
    fields,
    center="c",
    width=None,
    axes_unit=None,
    origin="center-window",
    fontsize=18,
    field_parameters=None,
    window_size=8.0,
    aspect=None,
    data_source=None,
):
    r"""Creates a plot of a 2D dataset

    Given a ds object and a field name string, this will return a
    PWViewerMPL object containing the plot.

    The plot can be updated using one of the many helper functions
    defined in PlotWindow.

    Parameters
    ----------
    ds : `Dataset`
         This is the dataset object corresponding to the
         simulation output to be plotted.
    fields : string
         The name of the field(s) to be plotted.
    center : A sequence of floats, a string, or a tuple.
         The coordinate of the center of the image. If set to 'c', 'center' or
         left blank, the plot is centered on the middle of the domain. If set to
         'max' or 'm', the center will be located at the maximum of the
         ('gas', 'density') field. Centering on the max or min of a specific
         field is supported by providing a tuple such as ("min","temperature") or
         ("max","dark_matter_density"). Units can be specified by passing in *center*
         as a tuple containing a coordinate and string unit name or by passing
         in a YTArray. If a list or unitless array is supplied, code units are
         assumed. For plot_2d, this keyword accepts a coordinate in two dimensions.
    width : tuple or a float.
         Width can have four different formats to support windows with variable
         x and y widths.  They are:

         ==================================     =======================
         format                                 example
         ==================================     =======================
         (float, string)                        (10,'kpc')
         ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
         float                                  0.2
         (float, float)                         (0.2, 0.3)
         ==================================     =======================

         For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
         wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
         window that is 10 kiloparsecs wide along the x axis and 15
         kiloparsecs wide along the y axis.  In the other two examples, code
         units are assumed, for example (0.2, 0.3) requests a plot that has an
         x width of 0.2 and a y width of 0.3 in code units.  If units are
         provided the resulting plot axis labels will use the supplied units.
    origin : string or length 1, 2, or 3 sequence.
         The location of the origin of the plot coordinate system. This
         is typically represented by a '-' separated string or a tuple of
         strings. In the first index the y-location is given by 'lower',
         'upper', or 'center'. The second index is the x-location, given as
         'left', 'right', or 'center'. Finally, whether the origin is
         applied in 'domain' space, plot 'window' space or 'native'
         simulation coordinate system is given. For example, both
         'upper-right-domain' and ['upper', 'right', 'domain'] place the
         origin in the upper right hand corner of domain space. If x or y
         are not given, a value is inferred. For instance, 'left-domain'
         corresponds to the lower-left hand corner of the simulation domain,
         'center-domain' corresponds to the center of the simulation domain,
         or 'center-window' for the center of the plot window. In the event
         that none of these options place the origin in a desired location,
         a sequence of tuples and a string specifying the
         coordinate space can be given. If plain numeric types are input,
         units of `code_length` are assumed. Further examples:

         =============================================== ===============================
         format                                          example
         =============================================== ===============================
         '{space}'                                       'domain'
         '{xloc}-{space}'                                'left-window'
         '{yloc}-{space}'                                'upper-domain'
         '{yloc}-{xloc}-{space}'                         'lower-right-window'
         ('{space}',)                                    ('window',)
         ('{xloc}', '{space}')                           ('right', 'domain')
         ('{yloc}', '{space}')                           ('lower', 'window')
         ('{yloc}', '{xloc}', '{space}')                 ('lower', 'right', 'window')
         ((yloc, '{unit}'), (xloc, '{unit}'), '{space}') ((0, 'm'), (.4, 'm'), 'window')
         (xloc, yloc, '{space}')                         (0.23, 0.5, 'domain')
         =============================================== ===============================
    axes_unit : string
         The name of the unit for the tick labels on the x and y axes.
         Defaults to None, which automatically picks an appropriate unit.
         If axes_unit is '1', 'u', or 'unitary', it will not display the
         units, and only show the axes name.
    fontsize : integer
         The size of the fonts for the axis, colorbar, and tick labels.
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source: YTSelectionContainer object
         Object to be used for data selection. Defaults to ds.all_data(), a
         region covering the full domain
    """
    if ds.dimensionality != 2:
        raise RuntimeError("plot_2d only plots 2D datasets!")
    if ds.geometry in ["cartesian", "polar", "spectral_cube"]:
        axis = "z"
    elif ds.geometry == "cylindrical":
        axis = "theta"
    elif ds.geometry == "spherical":
        axis = "phi"
    else:
        raise NotImplementedError(
            f"plot_2d does not yet support datasets with {ds.geometry} geometries"
        )
    # Part of the convenience of plot_2d is to eliminate the use of the
    # superfluous coordinate, so we do that also with the center argument
    if not isinstance(center, str) and obj_length(center) == 2:
        c0_string = isinstance(center[0], str)
        c1_string = isinstance(center[1], str)
        if not c0_string and not c1_string:
            if obj_length(center[0]) == 2 and c1_string:
                center = ds.arr(center[0], center[1])
            elif not isinstance(center, YTArray):
                center = ds.arr(center, "code_length")
            center.convert_to_units("code_length")
        center = ds.arr([center[0], center[1], ds.domain_center[2]])
    return AxisAlignedSlicePlot(
        ds,
        axis,
        fields,
        center=center,
        width=width,
        axes_unit=axes_unit,
        origin=origin,
        fontsize=fontsize,
        field_parameters=field_parameters,
        window_size=window_size,
        aspect=aspect,
        data_source=data_source,
    )


def _check_right_handed(right_handed: bool, flip_horizontal: bool) -> bool:
    # temporary function to check if right_handed kwarg has been set. can
    # remove this after full depreciation.
    if not right_handed:
        issue_deprecation_warning(
            "The 'right_handed' argument is deprecated. Use 'flip_horizontal' callback instead.",
            since="4.1.0",
            removal="4.2.0",
        )
        if flip_horizontal:
            # may be able to remove this now that its not a kwarg
            raise ValueError(
                "Cannot use both 'right_handed' and 'flip_horizontal' arguments"
            )
    return right_handed
