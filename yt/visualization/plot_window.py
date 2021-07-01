from collections import defaultdict
from functools import wraps
from numbers import Number

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from more_itertools import always_iterable
from mpl_toolkits.axes_grid1 import ImageGrid
from packaging.version import parse as parse_version
from unyt.exceptions import UnitConversionError

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.config import ytcfg
from yt.data_objects.image_array import ImageArray
from yt.frontends.ytdata.data_structures import YTSpatialPlotDataset
from yt.funcs import fix_axis, fix_unitary, is_sequence, iter_fields, mylog, obj_length
from yt.units.unit_object import Unit
from yt.units.unit_registry import UnitParseError
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

from .base_plot_types import CallbackWrapper, ImagePlotMPL
from .fixed_resolution import (
    FixedResolutionBuffer,
    OffAxisProjectionFixedResolutionBuffer,
)
from .geo_plot_utils import get_mpl_transform
from .plot_container import (
    ImagePlotContainer,
    apply_callback,
    get_log_minorticks,
    get_symlog_minorticks,
    invalidate_data,
    invalidate_figure,
    invalidate_plot,
    linear_transform,
    log_transform,
    symlog_transform,
)
from .plot_modifications import callback_registry

import sys  # isort: skip

if sys.version_info < (3, 10):
    # this function is deprecated in more_itertools
    # because it is superseded by the standard library
    from more_itertools import zip_equal
else:

    def zip_equal(*args):
        # FUTURE: when only Python 3.10+ is supported,
        # drop this conditional and call the builtin zip
        # function directly where due
        return zip(*args, strict=True)


MPL_VERSION = parse_version(matplotlib.__version__)

# Some magic for dealing with pyparsing being included or not
# included in matplotlib (not in gentoo, yes in everything else)
try:
    from matplotlib.pyparsing_py3 import ParseFatalException
except ImportError:
    from pyparsing import ParseFatalException


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


class PlotWindow(ImagePlotContainer):
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
    ):
        self.center = None
        self._periodic = periodic
        self.oblique = oblique
        self._right_handed = right_handed
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

        for field in self.data_source._determine_fields(self.fields):
            finfo = self.data_source.ds._get_field_info(*field)
            if finfo.take_log:
                self._field_transform[field] = log_transform
            else:
                self._field_transform[field] = linear_transform

            log, linthresh = self._log_config[field]
            if log is not None:
                self.set_log(field, log, linthresh=linthresh)

            # Access the dictionary to force the key to be created
            self._units_config[field]

        self.setup_callbacks()
        self._setup_plots()

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

    def frb():
        doc = "The frb property."

        def fget(self):
            if self._frb is None or not self._data_valid:
                self._recreate_frb()
            return self._frb

        def fset(self, value):
            self._frb = value
            self._data_valid = True

        def fdel(self):
            del self._frb
            self._frb = None
            self._data_valid = False

        return locals()

    frb = property(**frb())

    def _recreate_frb(self):
        old_fields = None
        # If we are regenerating an frb, we want to know what fields we had before
        if self._frb is not None:
            old_fields = list(self._frb.keys())
            old_units = [str(self._frb[of].units) for of in old_fields]

        # Set the bounds
        if hasattr(self, "zlim"):
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
        )

        # At this point the frb has the valid bounds, size, aliasing, etc.
        if old_fields is None:
            self._frb._get_data_source_fields()

            # New frb, apply default units (if any)
            for field, field_unit in self._units_config.items():
                if field_unit is None:
                    continue

                field_unit = Unit(field_unit, registry=self.ds.unit_registry)
                is_projected = getattr(self, "projected", False)
                if is_projected:
                    # Obtain config
                    path_length_units = Unit(
                        ytcfg.get_most_specific(
                            "plot", *field, "path_length_units", fallback="cm"
                        ),
                        registry=self.ds.unit_registry,
                    )
                    units = field_unit * path_length_units
                else:
                    units = field_unit
                try:
                    self.frb[field].convert_to_units(units)
                except UnitConversionError:
                    msg = (
                        "Could not apply default units from configuration.\n"
                        "Tried converting projected field %s from %s to %s, retaining units %s:\n"
                        "\tgot units for field: %s"
                    )
                    args = [
                        field,
                        self.frb[field].units,
                        units,
                        field_unit,
                        units,
                    ]
                    if is_projected:
                        msg += "\n\tgot units for integration length: %s"
                        args += [path_length_units]

                    msg += "\nCheck your configuration file."

                    mylog.error(msg, *args)
        else:
            # Restore the old fields
            for key, units in zip(old_fields, old_units):
                self._frb[key]
                equiv = self._equivalencies[key]
                if equiv[0] is None:
                    self._frb[key].convert_to_units(units)
                else:
                    self.frb.set_unit(key, units, equiv[0], equiv[1])

        # Restore the override fields
        for key in self.override_fields:
            self._frb[key]

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
            raise RuntimeError(
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
            raise RuntimeError(
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
           The name of the new unit.

        equivalency : string, optional
           If set, the equivalency to use to convert the current units to
           the new requested unit. If None, the unit conversion will be done
           without an equivalency

        equivalency_kwargs : string, optional
           Keyword arguments to be passed to the equivalency. Only used if
           ``equivalency`` is set.
        """
        if equivalency_kwargs is None:
            equivalency_kwargs = {}
        field = self.data_source._determine_fields(field)[0]
        for f, u in zip_equal(iter_fields(field), always_iterable(new_unit)):
            self.frb.set_unit(f, u, equivalency, equivalency_kwargs)
            self._equivalencies[f] = (equivalency, equivalency_kwargs)
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

    def set_window_size(self, size):
        """This calls set_figure_size to adjust the size of the plot window."""
        from yt._maintenance.deprecation import issue_deprecation_warning

        issue_deprecation_warning(
            "`PlotWindow.set_window_size` is a deprecated alias "
            "for `PlotWindow.set_figure_size`.",
            removal="4.1.0",
        )
        self.set_figure_size(size)
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
        self._right_handed = not self._right_handed

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
    _frb_generator = None
    _plot_type = None
    _data_valid = False

    def __init__(self, *args, **kwargs):
        if self._frb_generator is None:
            self._frb_generator = kwargs.pop("frb_generator")
        if self._plot_type is None:
            self._plot_type = kwargs.pop("plot_type")
        self._splat_color = kwargs.pop("splat_color", None)
        PlotWindow.__init__(self, *args, **kwargs)

    def _setup_origin(self):
        origin = self.origin
        axis_index = self.data_source.axis
        xc = None
        yc = None

        if isinstance(origin, str):
            origin = tuple(origin.split("-"))[:3]
        if 1 == len(origin):
            origin = ("lower", "left") + origin
        elif 2 == len(origin) and origin[0] in {"left", "right", "center"}:
            o0map = {"left": "lower", "right": "upper", "center": "center"}
            origin = (o0map[origin[0]],) + origin
        elif 2 == len(origin) and origin[0] in {"lower", "upper", "center"}:
            origin = (origin[0], "center", origin[-1])
        elif 3 == len(origin) and isinstance(origin[0], (int, float)):
            xc = self.ds.quan(origin[0], "code_length")
            yc = self.ds.quan(origin[1], "code_length")
        elif 3 == len(origin) and isinstance(origin[0], tuple):
            xc = self.ds.quan(origin[0][0], origin[0][1])
            yc = self.ds.quan(origin[1][0], origin[0][1])

        assert origin[-1] in ["window", "domain", "native"]

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
        else:
            mylog.warning("origin = %s", origin)
            msg = (
                'origin keyword "{}" not recognized, must declare "domain" '
                'or "center" as the last term in origin.'
            ).format(self.origin)
            raise RuntimeError(msg)
        if xc is None and yc is None:
            if origin[0] == "lower":
                yc = yllim
            elif origin[0] == "upper":
                yc = yrlim
            elif origin[0] == "center":
                yc = (yllim + yrlim) / 2.0
            else:
                mylog.warning("origin = %s", origin)
                msg = (
                    'origin keyword "{0}" not recognized, must declare "lower" '
                    '"upper" or "center" as the first term in origin.'
                )
                msg = msg.format(self.origin)
                raise RuntimeError(msg)

            if origin[1] == "left":
                xc = xllim
            elif origin[1] == "right":
                xc = xrlim
            elif origin[1] == "center":
                xc = (xllim + xrlim) / 2.0
            else:
                mylog.warning("origin = %s", origin)
                msg = (
                    'origin keyword "{0}" not recognized, must declare "left" '
                    '"right" or "center" as the second term in origin.'
                )
                msg = msg.format(self.origin)
                raise RuntimeError(msg)

        x_in_bounds = xc >= xllim and xc <= xrlim
        y_in_bounds = yc >= yllim and yc <= yrlim

        if not x_in_bounds and not y_in_bounds:
            msg = "origin inputs not in bounds of specified coordinate system domain."
            msg = msg.format(self.origin)
            raise RuntimeError(msg)

        return xc, yc

    def _setup_plots(self):
        from matplotlib.mathtext import MathTextParser

        if self._plot_valid:
            return
        if not self._data_valid:
            self._recreate_frb()
            self._data_valid = True
        self._colorbar_valid = True
        for f in list(set(self.data_source._determine_fields(self.fields))):
            axis_index = self.data_source.axis

            xc, yc = self._setup_origin()
            if self.ds.unit_system._code_flag or self.ds.no_cgs_equiv_length:
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
                (unit_x, unit_y) = (unit, unit)
            else:
                (unit_x, unit_y) = self._axes_unit_names

            # For some plots we may set aspect by hand, such as for spectral cube data.
            # This will likely be replaced at some point by the coordinate handler
            # setting plot aspect.
            if self.aspect is None:
                self.aspect = float(
                    (self.ds.quan(1.0, unit_y) / self.ds.quan(1.0, unit_x)).in_cgs()
                )
            extentx = [(self.xlim[i] - xc).in_units(unit_x) for i in (0, 1)]
            extenty = [(self.ylim[i] - yc).in_units(unit_y) for i in (0, 1)]

            extent = extentx + extenty

            if f in self.plots.keys():
                zlim = (self.plots[f].zmin, self.plots[f].zmax)
            else:
                zlim = (None, None)

            image = self.frb[f]
            if self._field_transform[f] == log_transform:
                msg = None
                use_symlog = False
                if zlim != (None, None):
                    pass
                elif np.nanmax(image) == np.nanmin(image):
                    msg = f"Plotting {f}: All values = {np.nanmax(image)}"
                elif np.nanmax(image) <= 0:
                    msg = (
                        f"Plotting {f}: All negative values. Max = {np.nanmax(image)}."
                    )
                    use_symlog = True
                elif not np.any(np.isfinite(image)):
                    msg = f"Plotting {f}: All values = NaN."
                elif np.nanmax(image) > 0.0 and np.nanmin(image) < 0:
                    msg = (
                        f"Plotting {f}: Both positive and negative values. "
                        f"Min = {np.nanmin(image)}, Max = {np.nanmax(image)}."
                    )
                    use_symlog = True
                elif np.nanmax(image) > 0.0 and np.nanmin(image) == 0:
                    # normally, a LogNorm scaling would still be OK here because
                    # LogNorm will mask 0 values when calculating vmin. But
                    # due to a bug in matplotlib's imshow, if the data range
                    # spans many orders of magnitude while containing zero points
                    # vmin can get rescaled to 0, resulting in an error when the image
                    # gets drawn. So here we switch to symlog to avoid that until
                    # a fix is in -- see PR #3161 and linked issue.
                    cutoff_sigdigs = 15
                    if (
                        np.log10(np.nanmax(image[np.isfinite(image)]))
                        - np.log10(np.nanmin(image[image > 0]))
                        > cutoff_sigdigs
                    ):
                        msg = f"Plotting {f}: Wide range and zeros."
                        use_symlog = True
                if msg is not None:
                    mylog.warning(msg)
                    if use_symlog:
                        mylog.warning("Switching to symlog colorbar scaling.")
                        self._field_transform[f] = symlog_transform
                        self._field_transform[f].func = None
                    else:
                        mylog.warning("Switching to linear colorbar scaling.")
                        self._field_transform[f] = linear_transform

            font_size = self._font_properties.get_size()

            fig = None
            axes = None
            cax = None
            draw_colorbar = True
            draw_axes = True
            draw_frame = draw_axes
            if f in self.plots:
                draw_colorbar = self.plots[f]._draw_colorbar
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
            self.plots[f] = WindowPlotMPL(
                ia,
                self._field_transform[f].name,
                self._field_transform[f].func,
                self._colormap_config[f],
                extent,
                zlim,
                self.figure_size,
                font_size,
                self.aspect,
                fig,
                axes,
                cax,
                self._projection,
                self._transform,
            )

            if not self._right_handed:
                ax = self.plots[f].axes
                ax.invert_xaxis()

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
                    self.plots[f].image.set_extent((xmin, xmax, ymin, ymax))
                    self.plots[f].axes.set_aspect("auto")

            x_label, y_label, colorbar_label = self._get_axes_labels(f)

            if x_label is not None:
                labels[0] = x_label
            if y_label is not None:
                labels[1] = y_label

            self.plots[f].axes.set_xlabel(labels[0])
            self.plots[f].axes.set_ylabel(labels[1])

            color = self._background_color[f]

            self.plots[f].axes.set_facecolor(color)

            # Determine the units of the data
            units = Unit(self.frb[f].units, registry=self.ds.unit_registry)
            units = units.latex_representation()

            if colorbar_label is None:
                colorbar_label = image.info["label"]
                if hasattr(self, "projected"):
                    colorbar_label = "$\\rm{Projected }$ %s" % colorbar_label
                if units is None or units == "":
                    pass
                else:
                    colorbar_label += r"$\ \ \left(" + units + r"\right)$"

            parser = MathTextParser("Agg")
            try:
                parser.parse(colorbar_label)
            except ParseFatalException as err:
                raise YTCannotParseUnitDisplayName(f, colorbar_label, str(err))

            self.plots[f].cb.set_label(colorbar_label)

            # x-y axes minorticks
            if f not in self._minorticks:
                self._minorticks[f] = True
            if self._minorticks[f]:
                self.plots[f].axes.minorticks_on()
            else:
                self.plots[f].axes.minorticks_off()

            # colorbar minorticks
            if f not in self._cbar_minorticks:
                self._cbar_minorticks[f] = True

            if self._cbar_minorticks[f]:
                vmin = np.float64(self.plots[f].cb.norm.vmin)
                vmax = np.float64(self.plots[f].cb.norm.vmax)

                if self._field_transform[f] == linear_transform:
                    self.plots[f].cax.minorticks_on()

                elif self._field_transform[f] == symlog_transform:
                    flinthresh = 10 ** np.floor(
                        np.log10(self.plots[f].cb.norm.linthresh)
                    )
                    mticks = self.plots[f].image.norm(
                        get_symlog_minorticks(flinthresh, vmin, vmax)
                    )
                    self.plots[f].cax.yaxis.set_ticks(mticks, minor=True)

                elif self._field_transform[f] == log_transform:
                    if MPL_VERSION >= parse_version("3.0.0"):
                        self.plots[f].cax.minorticks_on()
                        self.plots[f].cax.xaxis.set_visible(False)
                    else:
                        mticks = self.plots[f].image.norm(
                            get_log_minorticks(vmin, vmax)
                        )
                        self.plots[f].cax.yaxis.set_ticks(mticks, minor=True)

                else:
                    mylog.error(
                        "Unable to draw cbar minorticks for field "
                        "%s with transform %s ",
                        f,
                        self._field_transform[f],
                    )
                    self._cbar_minorticks[f] = False

            if not self._cbar_minorticks[f]:
                self.plots[f].cax.minorticks_off()

            if not draw_axes:
                self.plots[f]._toggle_axes(draw_axes, draw_frame)

            if not draw_colorbar:
                self.plots[f]._toggle_colorbar(draw_colorbar)

        self._set_font_properties()
        self.run_callbacks()
        self._plot_valid = True

    def setup_callbacks(self):
        for key in callback_registry:
            ignored = ["PlotCallback"]
            if self._plot_type.startswith("OffAxis"):
                ignored += [
                    "ParticleCallback",
                    "ClumpContourCallback",
                    "GridBoundaryCallback",
                ]
            if self._plot_type == "OffAxisProjection":
                ignored += [
                    "VelocityCallback",
                    "MagFieldCallback",
                    "QuiverCallback",
                    "CuttingQuiverCallback",
                    "StreamlineCallback",
                ]
            if self._plot_type == "Particle":
                ignored += [
                    "HopCirclesCallback",
                    "HopParticleCallback",
                    "ClumpContourCallback",
                    "GridBoundaryCallback",
                    "VelocityCallback",
                    "MagFieldCallback",
                    "QuiverCallback",
                    "CuttingQuiverCallback",
                    "StreamlineCallback",
                    "ContourCallback",
                ]
            if key in ignored:
                continue
            cbname = callback_registry[key]._type_name

            # We need to wrap to create a closure so that
            # CallbackMaker is bound to the wrapped method.
            def closure():
                CallbackMaker = callback_registry[key]

                @wraps(CallbackMaker)
                def method(*args, **kwargs):
                    # We need to also do it here as "invalidate_plot"
                    # and "apply_callback" require the functions'
                    # __name__ in order to work properly
                    @wraps(CallbackMaker)
                    def cb(self, *a, **kwa):
                        # We construct the callback method
                        # skipping self
                        return CallbackMaker(*a, **kwa)

                    # Create callback
                    cb = invalidate_plot(apply_callback(cb))

                    return cb(self, *args, **kwargs)

                return method

            self.__dict__["annotate_" + cbname] = closure()

    def annotate_clear(self, index=None):
        """
        Clear callbacks from the plot.  If index is not set, clear all
        callbacks.  If index is set, clear that index (ie 0 is the first one
        created, 1 is the 2nd one created, -1 is the last one created, etc.)

        .. note::

            Deprecated in favor of `clear_annotations`.

        See Also
        --------
        :py:meth:`yt.visualization.plot_window.PWViewerMPL.clear_annotations`
        """
        issue_deprecation_warning(
            "`annotate_clear` has been deprecated "
            "in favor of `clear_annotations`. Using `clear_annotations`.",
            since="4.0.0",
            removal="4.1.0",
        )
        self.clear_annotations(index=index)

    @invalidate_plot
    def clear_annotations(self, index=None):
        """
        Clear callbacks from the plot.  If index is not set, clear all
        callbacks.  If index is set, clear that index (ie 0 is the first one
        created, 1 is the 2nd one created, -1 is the last one created, etc.)
        """
        if index is None:
            self._callbacks = []
        else:
            del self._callbacks[index]
        self.setup_callbacks()
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
            for name, (args, kwargs) in self._callbacks:
                cbw = CallbackWrapper(
                    self,
                    self.plots[f],
                    self.frb,
                    f,
                    self._font_properties,
                    self._font_color,
                )
                CallbackMaker = callback_registry[name]
                callback = CallbackMaker(*args[1:], **kwargs)
                try:
                    callback(cbw)
                except YTDataTypeUnsupported as e:
                    raise e
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


class AxisAlignedSlicePlot(PWViewerMPL):
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
    axis : int or one of 'x', 'y', 'z'
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
         Whether the implicit east vector for the image generated is set to make a right
         handed coordinate system with a normal vector, the direction of the
         'window' into the data.
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
        axis,
        fields,
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
    ):
        # this will handle time series data and controllers
        axis = fix_axis(axis, ds)
        (bounds, center, display_center) = get_window_parameters(
            axis, center, width, ds
        )
        if field_parameters is None:
            field_parameters = {}

        if ds.geometry in (
            "spherical",
            "cylindrical",
            "geographic",
            "internal_geographic",
        ):
            mylog.info("Setting origin='native' for %s geometry.", ds.geometry)
            origin = "native"

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
        )
        if axes_unit is None:
            axes_unit = get_axes_unit(width, ds)
        self.set_axes_unit(axes_unit)


class ProjectionPlot(PWViewerMPL):
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
    axis : int or one of 'x', 'y', 'z'
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
         Whether the implicit east vector for the image generated is set to make a right
         handed coordinate system with the direction of the
         'window' into the data.
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

         "mip" : pick out the maximum value of the field in the line of sight.

         "sum" : This method is the same as integrate, except that it does not
         multiply by a path length when performing the integration, and is
         just a straight summation of the field along the given axis. WARNING:
         This should only be used for uniform resolution grid datasets, as other
         datasets may result in unphysical images.
    proj_style : string
         The method of projection--same as method keyword.  Deprecated as of
         version 3.0.2.  Please use method instead.
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
         used.  Effectively sets a resolution limit to the image if buff_size is
         smaller than the finest gridding.

    Examples
    --------

    Create a projection plot with a width of 20 kiloparsecs centered on the
    center of the simulation box:

    >>> from yt import load
    >>> ds = load("IsolateGalaxygalaxy0030/galaxy0030")
    >>> p = ProjectionPlot(ds, "z", ("gas", "density"), width=(20, "kpc"))

    """
    _plot_type = "Projection"
    _frb_generator = FixedResolutionBuffer

    def __init__(
        self,
        ds,
        axis,
        fields,
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
        proj_style=None,
        window_size=8.0,
        buff_size=(800, 800),
        aspect=None,
    ):
        axis = fix_axis(axis, ds)
        if ds.geometry in (
            "spherical",
            "cylindrical",
            "geographic",
            "internal_geographic",
        ):
            mylog.info("Setting origin='native' for %s geometry.", ds.geometry)
            origin = "native"
        if proj_style is not None:
            issue_deprecation_warning(
                "`proj_style` parameter is deprecated, use `method` instead.",
                removal="4.1.0",
            )
            method = proj_style
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
            )
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
        )
        if axes_unit is None:
            axes_unit = get_axes_unit(width, ds)
        self.set_axes_unit(axes_unit)


class OffAxisSlicePlot(PWViewerMPL):
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
         Whether the implicit east vector for the image generated is set to make a right
         handed coordinate system with the north vector and the normal, the direction of
         the 'window' into the data.
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
    ):
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
    _key_fields = []

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
    ):
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

    def _determine_fields(self, *args):
        return self.dd._determine_fields(*args)


class OffAxisProjectionPlot(PWViewerMPL):
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
    depth : A tuple or a float
         A tuple containing the depth to project through and the string
         key of the unit: (width, 'unit').  If set to a float, code units
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
         A vector defining the 'up' direction in the plot.  This
         option sets the orientation of the slicing plane.  If not
         set, an arbitrary grid-aligned north-vector is chosen.
    right_handed : boolean
         Whether the implicit east vector for the image generated is set to make a right
         handed coordinate system with the north vector and the normal, the direction of
         the 'window' into the data.
    fontsize : integer
         The size of the fonts for the axis, colorbar, and tick labels.
    method : string
         The method of projection.  Valid methods are:

         "integrate" with no weight_field specified : integrate the requested
         field along the line of sight.

         "integrate" with a weight_field specified : weight the requested
         field by the weighting field and integrate along the line of sight.

         "sum" : This method is the same as integrate, except that it does not
         multiply by a path length when performing the integration, and is
         just a straight summation of the field along the given axis. WARNING:
         This should only be used for uniform resolution grid datasets, as other
         datasets may result in unphysical images.
    data_source: YTSelectionContainer object
         Object to be used for data selection. Defaults to ds.all_data(), a
         region covering the full domain
    buff_size: length 2 sequence
         Size of the buffer to use for the image, i.e. the number of resolution elements
         used.  Effectively sets a resolution limit to the image if buff_size is
         smaller than the finest gridding.
    """
    _plot_type = "OffAxisProjection"
    _frb_generator = OffAxisProjectionFixedResolutionBuffer

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
        data_source=None,
        buff_size=(800, 800),
    ):
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
        )

        validate_mesh_fields(OffAxisProj, fields)

        if max_level is not None:
            OffAxisProj.dd.max_level = max_level

        # If a non-weighted, integral projection, assure field label
        # reflects that
        if weight_field is None and OffAxisProj.method == "integrate":
            self.projected = True

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
        cbname,
        cblinthresh,
        cmap,
        extent,
        zlim,
        figure_size,
        fontsize,
        aspect,
        figure,
        axes,
        cax,
        mpl_proj,
        mpl_transform,
    ):
        from matplotlib.ticker import ScalarFormatter

        self._draw_colorbar = True
        self._draw_axes = True
        self._draw_frame = True
        self._fontsize = fontsize
        self._figure_size = figure_size
        self._projection = mpl_proj
        self._transform = mpl_transform

        # Compute layout
        fontscale = float(fontsize) / 18.0
        if fontscale < 1.0:
            fontscale = np.sqrt(fontscale)

        if is_sequence(figure_size):
            fsize = figure_size[0]
        else:
            fsize = figure_size
        self._cb_size = 0.0375 * fsize
        self._ax_text_size = [1.2 * fontscale, 0.9 * fontscale]
        self._top_buff_size = 0.30 * fontscale
        self._aspect = ((extent[1] - extent[0]) / (extent[3] - extent[2])).in_cgs()
        self._unit_aspect = aspect

        size, axrect, caxrect = self._get_best_layout()

        super().__init__(size, axrect, caxrect, zlim, figure, axes, cax)

        self._init_image(data, cbname, cblinthresh, cmap, extent, aspect)

        # In matplotlib 2.1 and newer we'll be able to do this using
        # self.image.axes.ticklabel_format
        # See https://github.com/matplotlib/matplotlib/pull/6337
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-2, 3))
        self.image.axes.xaxis.set_major_formatter(formatter)
        self.image.axes.yaxis.set_major_formatter(formatter)
        if cbname == "linear":
            self.cb.formatter.set_scientific(True)
            try:
                self.cb.formatter.set_useMathText(True)
            except AttributeError:
                # this is only available in mpl > 2.1
                pass
            self.cb.formatter.set_powerlimits((-2, 3))
            self.cb.update_ticks()

    def _create_axes(self, axrect):
        self.axes = self.figure.add_axes(axrect, projection=self._projection)


def SlicePlot(ds, normal=None, fields=None, axis=None, *args, **kwargs):
    r"""
    A factory function for
    :class:`yt.visualization.plot_window.AxisAlignedSlicePlot`
    and :class:`yt.visualization.plot_window.OffAxisSlicePlot` objects.  This
    essentially allows for a single entry point to both types of slice plots,
    the distinction being determined by the specified normal vector to the
    slice.

    The returned plot object can be updated using one of the many helper
    functions defined in PlotWindow.

    Parameters
    ----------

    ds : :class:`yt.data_objects.static_output.Dataset`
        This is the dataset object corresponding to the
        simulation output to be plotted.
    normal : int or one of 'x', 'y', 'z', or sequence of floats
        This specifies the normal vector to the slice.  If given as an integer
        or a coordinate string (0=x, 1=y, 2=z), this function will return an
        :class:`AxisAlignedSlicePlot` object.  If given as a sequence of floats,
        this is interpreted as an off-axis vector and an
        :class:`OffAxisSlicePlot` object is returned.
    fields : string
         The name of the field(s) to be plotted.
    axis : int or one of 'x', 'y', 'z'
         An int corresponding to the axis to slice along (0=x, 1=y, 2=z)
         or the axis name itself.  If specified, this will replace normal.


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

    Raises
    ------

    AssertionError
        If a proper normal axis is not specified via the normal or axis
        keywords, and/or if a field to plot is not specified.

    Examples
    --------

    >>> from yt import load
    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> slc = SlicePlot(ds, "x", ("gas", "density"), center=[0.2, 0.3, 0.4])

    >>> slc = SlicePlot(
    ...     ds, [0.4, 0.2, -0.1], ("gas", "pressure"), north_vector=[0.2, -0.3, 0.1]
    ... )

    """
    if axis is not None:
        issue_deprecation_warning(
            "SlicePlot's argument 'axis' is a deprecated alias for 'normal', it "
            "will be removed in a future version of yt.",
            since="4.0.0",
            removal="4.1.0",
        )
        if normal is not None:
            raise TypeError(
                "SlicePlot() received incompatible arguments 'axis' and 'normal'"
            )
        normal = axis

    # to keep positional ordering we had to make 'normal' and 'fields' keywords
    if normal is None:
        raise TypeError("Missing argument in SlicePlot(): 'normal'")

    if fields is None:
        raise TypeError("Missing argument in SlicePlot(): 'fields'")

    # use an AxisAlignedSlicePlot where possible, e.g.:
    # maybe someone passed normal=[0,0,0.2] when they should have just used "z"
    if is_sequence(normal) and not isinstance(normal, str):
        if np.count_nonzero(normal) == 1:
            normal = ("x", "y", "z")[np.nonzero(normal)[0][0]]
        else:
            normal = np.array(normal, dtype="float64")
            np.divide(normal, np.dot(normal, normal), normal)

    # by now the normal should be properly set to get either a On/Off Axis plot
    if is_sequence(normal) and not isinstance(normal, str):
        # OffAxisSlicePlot has hardcoded origin; remove it if in kwargs
        if "origin" in kwargs:
            mylog.warning(
                "Ignoring 'origin' keyword as it is ill-defined for "
                "an OffAxisSlicePlot object."
            )
            del kwargs["origin"]

        return OffAxisSlicePlot(ds, normal, fields, *args, **kwargs)
    else:
        # north_vector not used in AxisAlignedSlicePlots; remove it if in kwargs
        if "north_vector" in kwargs:
            mylog.warning(
                "Ignoring 'north_vector' keyword as it is ill-defined for "
                "an AxisAlignedSlicePlot object."
            )
            del kwargs["north_vector"]

        return AxisAlignedSlicePlot(ds, normal, fields, *args, **kwargs)


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
