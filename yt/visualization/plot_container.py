import __builtin__
import base64
import numpy as np
import matplotlib
import os
import types
from functools import wraps
from matplotlib.font_manager import FontProperties

from .tick_locators import LogLocator, LinearLocator
from .color_maps import yt_colormaps, is_colormap
from .plot_modifications import \
    callback_registry
from .plot_window import \
    CallbackWrapper
from .base_plot_types import CallbackWrapper

from yt.funcs import \
    defaultdict, get_image_suffix, get_ipython_api_version
from yt.utilities.definitions import axis_names
from yt.utilities.exceptions import \
    YTNotInsideNotebook


def invalidate_data(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        rv = f(*args, **kwargs)
        args[0]._data_valid = False
        args[0]._plot_valid = False
        if hasattr(args[0], '_recreate_frb'):
            args[0]._recreate_frb()
        if args[0]._initfinished:
            args[0]._setup_plots()
        return rv
    return newfunc


def invalidate_figure(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        rv = f(*args, **kwargs)
        for field in args[0].fields:
            args[0].plots[field].figure = None
            args[0].plots[field].axes = None
            args[0].plots[field].cax = None
        return rv
    return newfunc


def invalidate_plot(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        rv = f(*args, **kwargs)
        args[0]._plot_valid = False
        args[0]._setup_plots()
        return rv
    return newfunc


def apply_callback(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        #rv = f(*args[1:], **kwargs)
        args[0]._callbacks.append((f.__name__, (args, kwargs)))
        return args[0]
    return newfunc

field_transforms = {}


class FieldTransform(object):
    def __init__(self, name, func, locator):
        self.name = name
        self.func = func
        self.locator = locator
        field_transforms[name] = self

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

    def ticks(self, mi, ma):
        try:
            ticks = self.locator(mi, ma)
        except:
            ticks = []
        return ticks

log_transform = FieldTransform('log10', np.log10, LogLocator())
linear_transform = FieldTransform('linear', lambda x: x, LinearLocator())


class ImagePlotContainer(object):
    """A countainer for plots with colorbars.

    """
    _plot_type = None
    _plot_valid = False
    _colorbar_valid = False

    def __init__(self, data_source, fields, figure_size, fontsize):
        self.data_source = data_source
        self.figure_size = figure_size
        self.plots = {}
        self._callbacks = []
        self._field_transform = {}
        self._colormaps = defaultdict(lambda: 'algae')
        font_path = matplotlib.get_data_path() + '/fonts/ttf/STIXGeneral.ttf'
        self._font_properties = FontProperties(size=fontsize, fname=font_path)
        self._font_color = None

    @invalidate_plot
    def set_log(self, field, log):
        """set a field to log or linear.

        Parameters
        ----------
        field : string
            the field to set a transform
        log : boolean
            Log on/off.

        """
        if field == 'all':
            fields = self.plots.keys()
        else:
            fields = [field]
        for field in fields:
            if log:
                self._field_transform[field] = log_transform
            else:
                self._field_transform[field] = linear_transform
        return self

    def get_log(self, field):
        """get the transform type of a field.

        Parameters
        ----------
        field : string
            the field to get a transform

        """
        log = {}
        if field == 'all':
            fields = self.plots.keys()
        else:
            fields = [field]
        for field in fields:
            if self._field_transform[field] == log_transform:
                log[field] = True
            else:
                log[field] = False
        return log

    @invalidate_plot
    def set_transform(self, field, name):
        if name not in field_transforms:
            raise KeyError(name)
        self._field_transform[field] = field_transforms[name]
        return self

    @invalidate_plot
    def set_cmap(self, field, cmap_name):
        """set the colormap for one of the fields

        Parameters
        ----------
        field : string
            the field to set the colormap
            if field == 'all', applies to all plots.
        cmap_name : string
            name of the colormap

        """

        if field is 'all':
            fields = self.plots.keys()
        else:
            fields = [field]
        for field in fields:
            self._colorbar_valid = False
            self._colormaps[field] = cmap_name
        return self

    def set_zlim(self, field, zmin, zmax, dynamic_range=None):
        # Left blank to be overriden in subclasses
        pass

    def setup_callbacks(self):
        # Left blank to be overriden in subclasses
        pass

    def _setup_plots(self):
        # Left blank to be overriden in subclasses
        pass

    def _switch_pf(self, new_pf):
        ds = self.data_source
        name = ds._type_name
        kwargs = dict((n, getattr(ds, n)) for n in ds._con_args)
        new_ds = getattr(new_pf.h, name)(**kwargs)
        self.pf = new_pf
        self.data_source = new_ds
        self._data_valid = self._plot_valid = False
        self._recreate_frb()
        self._setup_plots()

    def __getitem__(self, item):
        return self.plots[item]

    @property
    def fields(self):
        return self._frb.data.keys() + self.override_fields

    def run_callbacks(self, f):
        keys = self._frb.keys()
        for name, (args, kwargs) in self._callbacks:
            cbw = CallbackWrapper(self, self.plots[f], self._frb, f)
            CallbackMaker = callback_registry[name]
            callback = CallbackMaker(*args[1:], **kwargs)
            callback(cbw)
        for key in self._frb.keys():
            if key not in keys:
                del self._frb[key]

    @invalidate_plot
    @invalidate_figure
    def set_font(self, font_dict=None):
        """set the font and font properties

        Parameters
        ----------
        font_dict : dict
        A dict of keyword parameters to be passed to
        :py:class:`matplotlib.font_manager.FontProperties`.

        Possible keys include
        * family - The font family. Can be serif, sans-serif, cursive, 'fantasy' or
          'monospace'.
        * style - The font style. Either normal, italic or oblique.
        * color - A valid color string like 'r', 'g', 'red', 'cobalt', and
          'orange'.
        * variant: Either normal or small-caps.
        * size: Either an relative value of xx-small, x-small, small, medium,
          large, x-large, xx-large or an absolute font size, e.g. 12
        * stretch: A numeric value in the range 0-1000 or one of
          ultra-condensed, extra-condensed, condensed, semi-condensed, normal,
          semi-expanded, expanded, extra-expanded or ultra-expanded
        * weight: A numeric value in the range 0-1000 or one of ultralight,
          light, normal, regular, book, medium, roman, semibold, demibold, demi,
          bold, heavy, extra bold, or black

        See the matplotlib font manager API documentation for more details.
        http://matplotlib.org/api/font_manager_api.html

        Notes
        -----
        Mathtext axis labels will only obey the `size` and `color` keyword.

        Examples
        --------
        This sets the font to be 24-pt, blue, sans-serif, italic, and
        bold-face.

        >>> slc = SlicePlot(pf, 'x', 'Density')
        >>> slc.set_font({'family':'sans-serif', 'style':'italic',
                          'weight':'bold', 'size':24, 'color':'blue'})

        """
        if font_dict is None:
            font_dict = {}
        if 'color' in font_dict:
            self._font_color = font_dict.pop('color')
        self._font_properties = \
            FontProperties(**font_dict)
        return self

    def set_cmap(self, field, cmap):
        """set the colormap for one of the fields

        Parameters
        ----------
        field : string
            the field to set a transform
            if field == 'all', applies to all plots.
        cmap : string
            name of the colormap

        """
        if field == 'all':
            fields = self.plots.keys()
        else:
            fields = [field]

        for field in fields:
            self._colorbar_valid = False
            self._colormaps[field] = cmap
            if isinstance(cmap, types.StringTypes):
                if str(cmap) in yt_colormaps:
                    cmap = yt_colormaps[str(cmap)]
                elif hasattr(matplotlib.cm, cmap):
                    cmap = getattr(matplotlib.cm, cmap)
            if not is_colormap(cmap) and cmap is not None:
                raise RuntimeError("Colormap '%s' does not exist!" % str(cmap))
            self.plots[field].image.set_cmap(cmap)
        return self

    @invalidate_plot
    @invalidate_figure
    def set_figure_size(self, size):
        """Sets a new figure size for the plot

        parameters
        ----------
        size : float
            The size of the figure on the longest axis (in units of inches),
            including the margins but not the colorbar.
        """
        self.figure_size = size
        return self

    def save(self, name=None, mpl_kwargs=None):
        """saves the plot to disk.

        Parameters
        ----------
        name : string
           The base of the filename.  If name is a directory or if name is not
           set, the filename of the dataset is used.
        mpl_kwargs : dict
           A dict of keyword arguments to be passed to matplotlib.

        >>> slc.save(mpl_kwargs={'bbox_inches':'tight'})

        """
        names = []
        if mpl_kwargs is None: mpl_kwargs = {}
        if name is None:
            name = str(self.pf)
        name = os.path.expanduser(name)
        if name[-1] == os.sep and not os.path.isdir(name):
            os.mkdir(name)
        if os.path.isdir(name):
            name = name + (os.sep if name[-1] != os.sep else '') + str(self.pf)
        suffix = get_image_suffix(name)
        if suffix != '':
            for k, v in self.plots.iteritems():
                names.append(v.save(name, mpl_kwargs))
            return names
        axis = axis_names[self.data_source.axis]
        weight = None
        type = self._plot_type
        if type in ['Projection', 'OffAxisProjection']:
            weight = self.data_source.weight_field
            if weight is not None:
                weight = weight.replace(' ', '_')
        if 'Cutting' in self.data_source.__class__.__name__:
            type = 'OffAxisSlice'
        for k, v in self.plots.iteritems():
            if axis:
                n = "%s_%s_%s_%s" % (name, type, axis, k.replace(' ', '_'))
            else:
                # for cutting planes
                n = "%s_%s_%s" % (name, type, k.replace(' ', '_'))
            if weight:
                n += "_%s" % (weight)
            names.append(v.save(n, mpl_kwargs))
        return names

    @invalidate_data
    def refresh(self):
        # invalidate_data will take care of everything
        return self

    def _send_zmq(self):
        try:
            # pre-IPython v1.0
            from IPython.zmq.pylab.backend_inline import send_figure as display
        except ImportError:
            # IPython v1.0+
            from IPython.core.display import display
        for k, v in sorted(self.plots.iteritems()):
            # Due to a quirk in the matplotlib API, we need to create
            # a dummy canvas variable here that is never used.
            canvas = FigureCanvasAgg(v.figure)  # NOQA
            display(v.figure)

    def show(self):
        r"""This will send any existing plots to the IPython notebook.
        function name.

        If yt is being run from within an IPython session, and it is able to
        determine this, this function will send any existing plots to the
        notebook for display.

        If yt can't determine if it's inside an IPython session, it will raise
        YTNotInsideNotebook.

        Examples
        --------

        >>> from yt.mods import SlicePlot
        >>> slc = SlicePlot(pf, "x", ["Density", "VelocityMagnitude"])
        >>> slc.show()

        """
        if "__IPYTHON__" in dir(__builtin__):
            api_version = get_ipython_api_version()
            if api_version in ('0.10', '0.11'):
                self._send_zmq()
            else:
                from IPython.display import display
                display(self)
        else:
            raise YTNotInsideNotebook

    def display(self, name=None, mpl_kwargs=None):
        """Will attempt to show the plot in in an IPython notebook.
        Failing that, the plot will be saved to disk."""
        try:
            return self.show()
        except YTNotInsideNotebook:
            return self.save(name=name, mpl_kwargs=mpl_kwargs)

    def _repr_html_(self):
        """Return an html representation of the plot object. Will display as a
        png for each WindowPlotMPL instance in self.plots"""
        ret = ''
        for field in self.plots:
            img = base64.b64encode(self.plots[field]._repr_png_())
            ret += '<img src="data:image/png;base64,%s"><br>' % img
        return ret
