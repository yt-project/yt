import abc
import base64
import builtins
import os
import sys
import warnings
from collections import defaultdict
from functools import wraps
from typing import Any, Dict, List, Optional, Tuple, Type, Union

import matplotlib
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from matplotlib.font_manager import FontProperties
from unyt.dimensions import length

from yt._maintenance.deprecation import issue_deprecation_warning
from yt._typing import Quantity
from yt.config import ytcfg
from yt.data_objects.time_series import DatasetSeries
from yt.funcs import ensure_dir, is_sequence, iter_fields
from yt.units.unit_object import Unit  # type: ignore
from yt.utilities.definitions import formatted_length_unit_names
from yt.utilities.exceptions import YTConfigurationError, YTNotInsideNotebook
from yt.visualization._commons import get_default_from_config
from yt.visualization._handlers import ColorbarHandler, NormHandler
from yt.visualization.base_plot_types import PlotMPL

from ._commons import (
    DEFAULT_FONT_PROPERTIES,
    invalidate_data,
    invalidate_figure,
    invalidate_plot,
    validate_image_name,
    validate_plot,
)

if sys.version_info >= (3, 8):
    from typing import Final, Literal
else:
    from typing_extensions import Final, Literal

latex_prefixes = {
    "u": r"\mu",
}


def apply_callback(f):
    issue_deprecation_warning(
        "The apply_callback decorator is not used in yt any more and "
        "will be removed in a future version. "
        "Please do not use it.",
        since="4.1",
    )

    @wraps(f)
    def newfunc(*args, **kwargs):
        args[0]._callbacks.append((f.__name__, (args, kwargs)))
        return args[0]

    return newfunc


def accepts_all_fields(func):
    """
    Decorate a function whose second argument is <field> and deal with the special case
    field == 'all', looping over all fields already present in the PlotContainer object.

    """
    # This is to be applied to PlotContainer class methods with the following signature:
    #
    # f(self, field, *args, **kwargs) -> self
    @wraps(func)
    def newfunc(self, field, *args, **kwargs):
        if field == "all":
            field = self.plots.keys()
        for f in self.data_source._determine_fields(field):
            func(self, f, *args, **kwargs)
        return self

    return newfunc


# define a singleton sentinel to be used as default value distinct from None
class Unset:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = object.__new__(cls)
        return cls._instance


UNSET: Final = Unset()


class PlotDictionary(defaultdict):
    def __getitem__(self, item):
        return defaultdict.__getitem__(
            self, self.data_source._determine_fields(item)[0]
        )

    def __setitem__(self, item, value):
        return defaultdict.__setitem__(
            self, self.data_source._determine_fields(item)[0], value
        )

    def __contains__(self, item):
        return defaultdict.__contains__(
            self, self.data_source._determine_fields(item)[0]
        )

    def __init__(self, data_source, default_factory=None):
        self.data_source = data_source
        return defaultdict.__init__(self, default_factory)


class PlotContainer(abc.ABC):
    """A container for generic plots"""

    _plot_dict_type: Type[PlotDictionary] = PlotDictionary
    _plot_type: Optional[str] = None
    _plot_valid = False

    _default_figure_size = tuple(matplotlib.rcParams["figure.figsize"])
    _default_font_size = 14.0

    def __init__(self, data_source, figure_size=None, fontsize: Optional[float] = None):
        self.data_source = data_source
        self.ds = data_source.ds
        self.ts = self._initialize_dataset(self.ds)
        self.plots = self.__class__._plot_dict_type(data_source)

        self._set_figure_size(figure_size)

        if fontsize is None:
            fontsize = self.__class__._default_font_size
        if sys.version_info >= (3, 9):
            font_dict = DEFAULT_FONT_PROPERTIES | {"size": fontsize}
        else:
            font_dict = {**DEFAULT_FONT_PROPERTIES, "size": fontsize}  # type:ignore

        self._font_properties = FontProperties(**font_dict)
        self._font_color = None
        self._xlabel = None
        self._ylabel = None
        self._minorticks: Dict[Tuple[str, str], bool] = {}

    @accepts_all_fields
    @invalidate_plot
    def set_log(
        self,
        field,
        log: Optional[bool] = None,
        *,
        linthresh: Optional[Union[float, Quantity, Literal["auto"]]] = None,
        symlog_auto: Optional[bool] = None,  # deprecated
    ):
        """set a field to log, linear, or symlog.

        Symlog scaling is a combination of linear and log, where from 0 to a
        threshold value, it operates as linear, and then beyond that it operates as
        log.  Symlog can also work with negative values in log space as well as
        negative and positive values simultaneously and symmetrically.  If symlog
        scaling is desired, please set log=True and either set symlog_auto=True or
        select a value for linthresh.

        Parameters
        ----------
        field : string
            the field to set a transform
            if field == 'all', applies to all plots.
        log : boolean, optional
            set log to True for log scaling, False for linear scaling.
        linthresh : float, (float, str), unyt_quantity, or 'auto', optional
            when using symlog scaling, linthresh is the value at which scaling
            transitions from linear to logarithmic.  linthresh must be positive.
            Note: setting linthresh will automatically enable symlog scale

        Note that *log* and *linthresh* are mutually exclusive arguments
        """
        if log is None and linthresh is None and symlog_auto is None:
            raise TypeError("set_log requires log or linthresh be set")

        if symlog_auto is not None:
            issue_deprecation_warning(
                "the symlog_auto argument is deprecated. Use linthresh='auto' instead",
                since="4.1",
                stacklevel=5,
            )
            if symlog_auto is True:
                linthresh = "auto"
            elif symlog_auto is False:
                pass
            else:
                raise TypeError(
                    "Received invalid value for parameter symlog_auto. "
                    f"Expected a boolean, got {symlog_auto!r}"
                )

        if log is not None and linthresh is not None:
            # we do not raise an error here for backward compatibility
            warnings.warn(
                f"log={log} has no effect because linthresh specified. Using symlog.",
                stacklevel=4,
            )

        pnh = self.plots[field].norm_handler

        if linthresh is not None:
            if isinstance(linthresh, str):
                if linthresh == "auto":
                    pnh.norm_type = SymLogNorm
                else:
                    raise ValueError(
                        "Expected a number, a unyt_quantity, a (float, 'unit') tuple, or 'auto'. "
                        f"Got linthresh={linthresh!r}"
                    )
            else:
                # pnh takes care of switching to symlog when linthresh is set
                pnh.linthresh = linthresh
        elif log is True:
            pnh.norm_type = LogNorm
        elif log is False:
            pnh.norm_type = Normalize
        else:
            raise TypeError(
                f"Could not parse arguments log={log!r}, linthresh={linthresh!r}"
            )

        return self

    def get_log(self, field):
        """get the transform type of a field.

        Parameters
        ----------
        field : string
            the field to get a transform
            if field == 'all', applies to all plots.

        """
        # devnote : accepts_all_fields decorator is not applicable here because
        # the return variable isn't self
        issue_deprecation_warning(
            "The get_log method is not reliable and is deprecated. "
            "Please do not rely on it.",
            since="4.1",
        )
        log = {}
        if field == "all":
            fields = list(self.plots.keys())
        else:
            fields = field
        for field in self.data_source._determine_fields(fields):
            pnh = self.plots[field].norm_handler
            if pnh.norm is not None:
                log[field] = type(pnh.norm) is LogNorm
            elif pnh.norm_type is not None:
                log[field] = pnh.norm_type is LogNorm
            else:
                # the NormHandler object has no constraints yet
                # so we'll assume defaults
                log[field] = True
        return log

    @invalidate_plot
    def set_transform(self, field, name: str):
        field = self.data_source._determine_fields(field)[0]
        pnh = self.plots[field].norm_handler
        pnh.norm_type = {
            "linear": Normalize,
            "log10": LogNorm,
            "symlog": SymLogNorm,
        }[name]
        return self

    @accepts_all_fields
    @invalidate_plot
    def set_norm(self, field, norm: Normalize):
        r"""
        Set a custom ``matplotlib.colors.Normalize`` to plot *field*.

        Any constraints previously set with `set_log`, `set_zlim` will be
        dropped.

        Note that any float value attached to *norm* (e.g. vmin, vmax,
        vcenter ...) will be read in the current displayed units, which can be
        controlled with the `set_unit` method.

        Parameters
        ----------
        field : str or tuple[str, str]
            if field == 'all', applies to all plots.
        norm : matplotlib.colors.Normalize
            see https://matplotlib.org/stable/tutorials/colors/colormapnorms.html
        """
        pnh = self.plots[field].norm_handler
        pnh.norm = norm
        return self

    @accepts_all_fields
    @invalidate_plot
    def set_minorticks(self, field, state):
        """Turn minor ticks on or off in the current plot.

        Displaying minor ticks reduces performance; turn them off
        using set_minorticks('all', False) if drawing speed is a problem.

        Parameters
        ----------
        field : string
            the field to remove minorticks
            if field == 'all', applies to all plots.
        state : bool
            the state indicating 'on' (True) or 'off' (False)

        """
        self._minorticks[field] = state
        return self

    @abc.abstractmethod
    def _setup_plots(self):
        # Left blank to be overridden in subclasses
        pass

    def render(self) -> None:
        r"""Render plots.
        This operation is expensive and usually doesn't need to be requested explicitly.
        In most cases, yt handles rendering automatically and delays it as much as possible
        to avoid redundant calls on each plot modification (e.g. via `annotate_*` methods).

        However, valid use cases of this method include:
        - fine control of render (and clear) operations when yt plots are combined with plot
          customizations other than plot callbacks (`annotate_*`)
        - testing
        """
        # this is a pure alias to the historic `_setup_plots` method
        # which preserves backward compatibility for extension code
        self._setup_plots()

    def _initialize_dataset(self, ts):
        if not isinstance(ts, DatasetSeries):
            if not is_sequence(ts):
                ts = [ts]
            ts = DatasetSeries(ts)
        return ts

    @invalidate_data
    def _switch_ds(self, new_ds, data_source=None):
        old_object = self.data_source
        name = old_object._type_name
        kwargs = {n: getattr(old_object, n) for n in old_object._con_args}
        kwargs["center"] = getattr(old_object, "center", None)
        if data_source is not None:
            if name != "proj":
                raise RuntimeError(
                    "The data_source keyword argument "
                    "is only defined for projections."
                )
            kwargs["data_source"] = data_source

        self.ds = new_ds

        # A _hack_ for ParticleProjectionPlots
        if name == "Particle":
            from yt.visualization.particle_plots import (
                ParticleAxisAlignedDummyDataSource,
            )

            new_object = ParticleAxisAlignedDummyDataSource(ds=self.ds, **kwargs)
        else:
            new_object = getattr(new_ds, name)(**kwargs)

        self.data_source = new_object

        for d in "xyz":
            lim_name = d + "lim"
            if hasattr(self, lim_name):
                lim = getattr(self, lim_name)
                lim = tuple(new_ds.quan(l.value, str(l.units)) for l in lim)
                setattr(self, lim_name, lim)
        self.plots.data_source = new_object
        self._colorbar_label.data_source = new_object
        self._setup_plots()

    @validate_plot
    def __getitem__(self, item):
        return self.plots[item]

    def _set_font_properties(self):
        for f in self.plots:
            self.plots[f]._set_font_properties(self._font_properties, self._font_color)

    @invalidate_plot
    @invalidate_figure
    def set_font(self, font_dict=None):
        """

        Set the font and font properties.

        Parameters
        ----------

        font_dict : dict
            A dict of keyword parameters to be passed to
            :class:`matplotlib.font_manager.FontProperties`.

            Possible keys include:

            * family - The font family. Can be serif, sans-serif, cursive,
              'fantasy' or 'monospace'.
            * style - The font style. Either normal, italic or oblique.
            * color - A valid color string like 'r', 'g', 'red', 'cobalt',
              and 'orange'.
            * variant - Either normal or small-caps.
            * size - Either a relative value of xx-small, x-small, small,
              medium, large, x-large, xx-large or an absolute font size, e.g. 12
            * stretch - A numeric value in the range 0-1000 or one of
              ultra-condensed, extra-condensed, condensed, semi-condensed,
              normal, semi-expanded, expanded, extra-expanded or ultra-expanded
            * weight - A numeric value in the range 0-1000 or one of ultralight,
              light, normal, regular, book, medium, roman, semibold, demibold,
              demi, bold, heavy, extra bold, or black

            See the matplotlib font manager API documentation for more details.
            https://matplotlib.org/stable/api/font_manager_api.html

        Notes
        -----

        Mathtext axis labels will only obey the `size` and `color` keyword.

        Examples
        --------

        This sets the font to be 24-pt, blue, sans-serif, italic, and
        bold-face.

        >>> slc = SlicePlot(ds, "x", "Density")
        >>> slc.set_font(
        ...     {
        ...         "family": "sans-serif",
        ...         "style": "italic",
        ...         "weight": "bold",
        ...         "size": 24,
        ...         "color": "blue",
        ...     }
        ... )

        """

        if font_dict is None:
            font_dict = {}
        if "color" in font_dict:
            self._font_color = font_dict.pop("color")
        # Set default values if the user does not explicitly set them.
        # this prevents reverting to the matplotlib defaults.
        if sys.version_info >= (3, 9):
            font_dict = DEFAULT_FONT_PROPERTIES | font_dict
        else:
            font_dict = {**DEFAULT_FONT_PROPERTIES, **font_dict}
        self._font_properties = FontProperties(**font_dict)
        return self

    def set_font_size(self, size):
        """Set the size of the font used in the plot

        This sets the font size by calling the set_font function.  See set_font
        for more font customization options.

        Parameters
        ----------
        size : float
        The absolute size of the font in points (1 pt = 1/72 inch).

        """
        return self.set_font({"size": size})

    def _set_figure_size(self, size):
        if size is None:
            self.figure_size = self.__class__._default_figure_size
        elif is_sequence(size):
            if len(size) != 2:
                raise TypeError(f"Expected a single float or a pair, got {size}")
            self.figure_size = float(size[0]), float(size[1])
        else:
            self.figure_size = float(size)

    @invalidate_plot
    @invalidate_figure
    def set_figure_size(self, size):
        """Sets a new figure size for the plot

        parameters
        ----------
        size : float, a sequence of two floats, or None
            The size of the figure (in units of inches),  including the margins
            but not the colorbar. If a single float is passed, it's interpreted
            as the size along the long axis.
            Pass None to reset
        """
        self._set_figure_size(size)
        return self

    @validate_plot
    def save(
        self,
        name: Optional[Union[str, List[str], Tuple[str, ...]]] = None,
        suffix: Optional[str] = None,
        mpl_kwargs: Optional[Dict[str, Any]] = None,
    ):
        """saves the plot to disk.

        Parameters
        ----------
        name : string or tuple, optional
           The base of the filename. If name is a directory or if name is not
           set, the filename of the dataset is used. For a tuple, the
           resulting path will be given by joining the elements of the
           tuple
        suffix : string, optional
           Specify the image type by its suffix. If not specified, the output
           type will be inferred from the filename. Defaults to '.png'.
        mpl_kwargs : dict, optional
           A dict of keyword arguments to be passed to matplotlib.

        >>> slc.save(mpl_kwargs={"bbox_inches": "tight"})

        """
        names = []
        if mpl_kwargs is None:
            mpl_kwargs = {}
        elif "format" in mpl_kwargs:
            new_suffix = mpl_kwargs.pop("format")
            if new_suffix != suffix:
                warnings.warn(
                    f"Overriding suffix {suffix!r} with mpl_kwargs['format'] = {new_suffix!r}. "
                    "Use the `suffix` argument directly to suppress this warning."
                )
            suffix = new_suffix

        if name is None:
            name = str(self.ds)
        elif isinstance(name, (list, tuple)):
            if not all(isinstance(_, str) for _ in name):
                raise TypeError(
                    f"Expected a single str or an iterable of str, got {name!r}"
                )
            name = os.path.join(*name)

        name = os.path.expanduser(name)

        parent_dir, _, prefix1 = name.replace(os.sep, "/").rpartition("/")
        parent_dir = parent_dir.replace("/", os.sep)

        if parent_dir and not os.path.isdir(parent_dir):
            ensure_dir(parent_dir)

        if name.endswith(("/", os.path.sep)):
            name = os.path.join(name, str(self.ds))

        new_name = validate_image_name(name, suffix)
        if new_name == name:
            for v in self.plots.values():
                out_name = v.save(name, mpl_kwargs)
                names.append(out_name)
            return names

        name = new_name
        prefix, suffix = os.path.splitext(name)

        if hasattr(self.data_source, "axis"):
            axis = self.ds.coordinates.axis_name.get(self.data_source.axis, "")
        else:
            axis = None
        weight = None
        stddev = None
        plot_type = self._plot_type
        if plot_type in ["Projection", "OffAxisProjection"]:
            weight = self.data_source.weight_field
            if weight is not None:
                weight = weight[1].replace(" ", "_")
            if getattr(self.data_source, "moment", 1) == 2:
                stddev = "standard_deviation"
        if "Cutting" in self.data_source.__class__.__name__:
            plot_type = "OffAxisSlice"

        for k, v in self.plots.items():
            if isinstance(k, tuple):
                k = k[1]

            if plot_type is None:
                # implemented this check to make mypy happy, because we can't use str.join
                # with PlotContainer._plot_type = None
                raise TypeError(f"{self.__class__} is missing a _plot_type value (str)")

            name_elements = [prefix, plot_type]
            if axis:
                name_elements.append(axis)
            name_elements.append(k.replace(" ", "_"))
            if weight:
                name_elements.append(weight)
            if stddev:
                name_elements.append(stddev)
            name = "_".join(name_elements) + suffix
            names.append(v.save(name, mpl_kwargs))
        return names

    @invalidate_data
    def refresh(self):
        # invalidate_data will take care of everything
        return self

    @validate_plot
    def show(self):
        r"""This will send any existing plots to the IPython notebook.

        If yt is being run from within an IPython session, and it is able to
        determine this, this function will send any existing plots to the
        notebook for display.

        If yt can't determine if it's inside an IPython session, it will raise
        YTNotInsideNotebook.

        Examples
        --------

        >>> from yt.mods import SlicePlot
        >>> slc = SlicePlot(
        ...     ds, "x", [("gas", "density"), ("gas", "velocity_magnitude")]
        ... )
        >>> slc.show()

        """
        interactivity = self.plots[list(self.plots.keys())[0]].interactivity
        if interactivity:
            for v in sorted(self.plots.values()):
                v.show()
        else:
            if "__IPYTHON__" in dir(builtins):
                from IPython.display import display

                display(self)
            else:
                raise YTNotInsideNotebook

    @validate_plot
    def display(self, name=None, mpl_kwargs=None):
        """Will attempt to show the plot in in an IPython notebook.
        Failing that, the plot will be saved to disk."""
        try:
            return self.show()
        except YTNotInsideNotebook:
            return self.save(name=name, mpl_kwargs=mpl_kwargs)

    @validate_plot
    def _repr_html_(self):
        """Return an html representation of the plot object. Will display as a
        png for each WindowPlotMPL instance in self.plots"""
        ret = ""
        for field in self.plots:
            img = base64.b64encode(self.plots[field]._repr_png_()).decode()
            ret += (
                r'<img style="max-width:100%;max-height:100%;" '
                r'src="data:image/png;base64,{}"><br>'.format(img)
            )
        return ret

    @invalidate_plot
    def set_xlabel(self, label):
        r"""
        Allow the user to modify the X-axis title
        Defaults to the global value. Fontsize defaults
        to 18.

        Parameters
        ----------
        label : str
            The new string for the x-axis.

        >>> plot.set_xlabel("H2I Number Density (cm$^{-3}$)")

        """
        self._xlabel = label
        return self

    @invalidate_plot
    def set_ylabel(self, label):
        r"""
        Allow the user to modify the Y-axis title
        Defaults to the global value.

        Parameters
        ----------
        label : str
            The new string for the y-axis.

        >>> plot.set_ylabel("Temperature (K)")

        """
        self._ylabel = label
        return self

    def _get_axes_unit_labels(self, unit_x, unit_y):
        axes_unit_labels = ["", ""]
        comoving = False
        hinv = False
        for i, un in enumerate((unit_x, unit_y)):
            unn = None
            if hasattr(self.data_source, "axis"):
                if hasattr(self.ds.coordinates, "image_units"):
                    # This *forces* an override
                    unn = self.ds.coordinates.image_units[self.data_source.axis][i]
                elif hasattr(self.ds.coordinates, "default_unit_label"):
                    axax = getattr(self.ds.coordinates, f"{'xy'[i]}_axis")[
                        self.data_source.axis
                    ]
                    unn = self.ds.coordinates.default_unit_label.get(axax, None)
            if unn in (1, "1", "dimensionless"):
                axes_unit_labels[i] = ""
                continue
            if unn is not None:
                axes_unit_labels[i] = r"\ \ \left(" + unn + r"\right)"
                continue
            # Use sympy to factor h out of the unit.  In this context 'un'
            # is a string, so we call the Unit constructor.
            expr = Unit(un, registry=self.ds.unit_registry).expr
            h_expr = Unit("h", registry=self.ds.unit_registry).expr
            # See http://docs.sympy.org/latest/modules/core.html#sympy.core.expr.Expr
            h_power = expr.as_coeff_exponent(h_expr)[1]
            # un is now the original unit, but with h factored out.
            un = str(expr * h_expr ** (-1 * h_power))
            un_unit = Unit(un, registry=self.ds.unit_registry)
            cm = Unit("cm").expr
            if str(un).endswith("cm") and cm not in un_unit.expr.atoms():
                comoving = True
                un = un[:-2]
            # no length units besides code_length end in h so this is safe
            if h_power == -1:
                hinv = True
            elif h_power != 0:
                # It doesn't make sense to scale a position by anything
                # other than h**-1
                raise RuntimeError
            if un not in ["1", "u", "unitary"]:
                if un in formatted_length_unit_names:
                    un = formatted_length_unit_names[un]
                else:
                    un = Unit(un, registry=self.ds.unit_registry)
                    un = un.latex_representation()
                    if hinv:
                        un = un + r"\,h^{-1}"
                    if comoving:
                        un = un + r"\,(1+z)^{-1}"
                    pp = un[0]
                    if pp in latex_prefixes:
                        symbol_wo_prefix = un[1:]
                        if symbol_wo_prefix in self.ds.unit_registry.prefixable_units:
                            un = un.replace(pp, "{" + latex_prefixes[pp] + "}", 1)
                if r"\frac" in un:
                    axes_unit_labels[i] = r"\ \ \left(" + un + r"\right)"
                else:
                    axes_unit_labels[i] = r"\ \ (" + un + r")"
        return axes_unit_labels

    def hide_colorbar(self, field=None):
        """
        Hides the colorbar for a plot and updates the size of the
        plot accordingly.  Defaults to operating on all fields for a
        PlotContainer object.

        Parameters
        ----------

        field : string, field tuple, or list of strings or field tuples (optional)
            The name of the field(s) that we want to hide the colorbar.
            If None or 'all' is provided, will default to using all fields available
            for this object.

        Examples
        --------

        This will save an image with no colorbar.

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> s = SlicePlot(ds, 2, "density", "c", (20, "kpc"))
        >>> s.hide_colorbar()
        >>> s.save()

        This will save an image with no axis or colorbar.

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> s = SlicePlot(ds, 2, "density", "c", (20, "kpc"))
        >>> s.hide_axes()
        >>> s.hide_colorbar()
        >>> s.save()
        """
        if field is None or field == "all":
            field = self.plots.keys()
        for f in self.data_source._determine_fields(field):
            self.plots[f].hide_colorbar()
        return self

    def show_colorbar(self, field=None):
        """
        Shows the colorbar for a plot and updates the size of the
        plot accordingly.  Defaults to operating on all fields for a
        PlotContainer object.  See hide_colorbar().

        Parameters
        ----------

        field : string, field tuple, or list of strings or field tuples (optional)
        The name of the field(s) that we want to show the colorbar.
        """
        if field is None:
            field = self.fields
        for f in iter_fields(field):
            self.plots[f].show_colorbar()
        return self

    def hide_axes(self, field=None, draw_frame=None):
        """
        Hides the axes for a plot and updates the size of the
        plot accordingly.  Defaults to operating on all fields for a
        PlotContainer object.

        Parameters
        ----------

        field : string, field tuple, or list of strings or field tuples (optional)
            The name of the field(s) that we want to hide the axes.

        draw_frame : boolean
            If True, the axes frame will still be drawn. Defaults to False.
            See note below for more details.

        Examples
        --------

        This will save an image with no axes.

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> s = SlicePlot(ds, 2, "density", "c", (20, "kpc"))
        >>> s.hide_axes()
        >>> s.save()

        This will save an image with no axis or colorbar.

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> s = SlicePlot(ds, 2, "density", "c", (20, "kpc"))
        >>> s.hide_axes()
        >>> s.hide_colorbar()
        >>> s.save()

        Note
        ----
        By default, when removing the axes, the patch on which the axes are
        drawn is disabled, making it impossible to later change e.g. the
        background colour. To force the axes patch to be displayed while still
        hiding the axes, set the ``draw_frame`` keyword argument to ``True``.
        """
        if field is None:
            field = self.fields
        for f in iter_fields(field):
            self.plots[f].hide_axes(draw_frame=draw_frame)
        return self

    def show_axes(self, field=None):
        """
        Shows the axes for a plot and updates the size of the
        plot accordingly.  Defaults to operating on all fields for a
        PlotContainer object.  See hide_axes().

        Parameters
        ----------

        field : string, field tuple, or list of strings or field tuples (optional)
            The name of the field(s) that we want to show the axes.
        """
        if field is None:
            field = self.fields
        for f in iter_fields(field):
            self.plots[f].show_axes()
        return self


class ImagePlotContainer(PlotContainer, abc.ABC):
    """A container for plots with colorbars."""

    _colorbar_valid = False

    def __init__(self, data_source, figure_size, fontsize):
        super().__init__(data_source, figure_size, fontsize)
        self._callbacks = []
        self._colorbar_label = PlotDictionary(self.data_source, lambda: None)

    def _get_default_handlers(
        self, field, default_display_units: Unit
    ) -> Tuple[NormHandler, ColorbarHandler]:

        usr_units_str = get_default_from_config(
            self.data_source, field=field, keys="units", defaults=[None]
        )
        if usr_units_str is not None:
            usr_units = Unit(usr_units_str)
            d1 = usr_units.dimensions
            d2 = default_display_units.dimensions

            if d1 == d2:
                display_units = usr_units
            elif getattr(self, "projected", False) and d2 / d1 == length:
                path_length_units = Unit(
                    ytcfg.get_most_specific(
                        "plot", *field, "path_length_units", fallback="cm"
                    ),
                    registry=self.data_source.ds.unit_registry,
                )
                display_units = usr_units * path_length_units
            else:
                raise YTConfigurationError(
                    f"Invalid units in configuration file for field {field!r}. "
                    f"Found {usr_units!r}"
                )
        else:
            display_units = default_display_units

        pnh = NormHandler(self.data_source, display_units=display_units)

        cbh = ColorbarHandler(
            cmap=get_default_from_config(
                self.data_source,
                field=field,
                keys="cmap",
                defaults=[None],
            )
        )
        return pnh, cbh

    @accepts_all_fields
    @invalidate_plot
    def set_cmap(self, field, cmap):
        """set the colormap for one of the fields

        Parameters
        ----------
        field : string
            the field to set the colormap
            if field == 'all', applies to all plots.
        cmap : string or tuple
            If a string, will be interpreted as name of the colormap.
            If a tuple, it is assumed to be of the form (name, type, number)
            to be used for palettable functionality. (name, type, number, bool)
            can be used to specify if a reverse colormap is to be used.

        """
        self._colorbar_valid = False
        self.plots[field].colorbar_handler.cmap = cmap
        return self

    @accepts_all_fields
    @invalidate_plot
    def set_background_color(self, field, color=None):
        """set the background color to match provided color

        Parameters
        ----------
        field : string
            the field to set the colormap
            if field == 'all', applies to all plots.
        color : string or RGBA tuple (optional)
            if set, set the background color to this color
            if unset, background color is set to the bottom value of
            the color map

        """
        cbh = self[field].colorbar_handler
        cbh.background_color = color
        return self

    @accepts_all_fields
    @invalidate_plot
    def set_zlim(
        self,
        field,
        zmin: Union[float, Quantity, Literal["min"], Unset] = UNSET,
        zmax: Union[float, Quantity, Literal["max"], Unset] = UNSET,
        dynamic_range: Optional[float] = None,
    ):
        """set the scale of the colormap

        Parameters
        ----------
        field : string
            the field to set a colormap scale
            if field == 'all', applies to all plots.
        zmin : float, Quantity, or 'min'
            the new minimum of the colormap scale. If 'min', will
            set to the minimum value in the current view.
        zmax : float, Quantity, or 'max'
            the new maximum of the colormap scale. If 'max', will
            set to the maximum value in the current view.

        Other Parameters
        ----------------
        dynamic_range : float (default: None)
            The dynamic range of the image.
            If zmin == None, will set zmin = zmax / dynamic_range
            If zmax == None, will set zmax = zmin * dynamic_range

        """
        if zmin is UNSET and zmax is UNSET:
            raise TypeError("Missing required argument zmin or zmax")

        if zmin is UNSET:
            zmin = None
        elif zmin is None:
            # this sentinel value juggling is barely maintainable
            # this use case is deprecated so we can simplify the logic here
            # in the future and use `None` as the default value,
            # instead of the custom sentinel UNSET
            issue_deprecation_warning(
                "Passing `zmin=None` explicitly is deprecated. "
                "If you wish to explicitly set zmin to the minimal "
                "data value, pass `zmin='min'` instead. "
                "Otherwise leave this argument unset.",
                since="4.1.0",
                stacklevel=5,
            )
            zmin = "min"

        if zmax is UNSET:
            zmax = None
        elif zmax is None:
            # see above
            issue_deprecation_warning(
                "Passing `zmax=None` explicitly is deprecated. "
                "If you wish to explicitly set zmax to the maximal "
                "data value, pass `zmin='max'` instead. "
                "Otherwise leave this argument unset.",
                since="4.1.0",
                stacklevel=5,
            )
            zmax = "max"

        pnh = self.plots[field].norm_handler
        pnh.vmin = zmin
        pnh.vmax = zmax
        pnh.dynamic_range = dynamic_range

        return self

    @accepts_all_fields
    @invalidate_plot
    def set_colorbar_minorticks(self, field, state):
        """turn colorbar minor ticks on or off in the current plot

        Displaying minor ticks reduces performance; turn them off
        using set_colorbar_minorticks('all', False) if drawing speed is a problem.

        Parameters
        ----------
        field : string
            the field to remove colorbar minorticks
            if field == 'all', applies to all plots.
        state : bool
            the state indicating 'on' (True) or 'off' (False)
        """
        self.plots[field].colorbar_handler.draw_minorticks = state
        return self

    @invalidate_plot
    def set_colorbar_label(self, field, label):
        r"""
        Sets the colorbar label.

        Parameters
        ----------
        field : str or tuple
          The name of the field to modify the label for.
        label : str
          The new label

        >>> plot.set_colorbar_label(
        ...     ("gas", "density"), "Dark Matter Density (g cm$^{-3}$)"
        ... )

        """
        field = self.data_source._determine_fields(field)
        self._colorbar_label[field] = label
        return self

    def _get_axes_labels(self, field):
        return (self._xlabel, self._ylabel, self._colorbar_label[field])


class BaseLinePlot(PlotContainer, abc.ABC):

    # A common ancestor to LinePlot and ProfilePlot

    @abc.abstractmethod
    def _get_axrect(self):
        pass

    def _get_plot_instance(self, field):
        if field in self.plots:
            return self.plots[field]
        axrect = self._get_axrect()

        pnh = NormHandler(self.data_source, display_units=self.data_source[field].units)
        finfo = self.data_source.ds._get_field_info(*field)
        if not finfo.take_log:
            pnh.norm_type = Normalize
        plot = PlotMPL(self.figure_size, axrect, norm_handler=pnh)
        self.plots[field] = plot

        return plot
