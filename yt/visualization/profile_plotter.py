"""
This is a simple mechanism for interfacing with Profile and Phase plots



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import base64
import types

from collections import defaultdict
from functools import wraps
from itertools import izip, repeat
import matplotlib
import numpy as np
import cStringIO

from matplotlib.font_manager import FontProperties

from .plot_window import WindowPlotMPL
from .base_plot_types import ImagePlotMPL
from .plot_container import ImagePlotContainer
from .image_writer import \
    write_image, apply_colormap
from yt.data_objects.profiles import \
     create_profile
from yt.utilities.lib import \
    write_png_to_string
from yt.data_objects.profiles import \
    BinnedProfile1D, \
    BinnedProfile2D
from yt.utilities.logger import ytLogger as mylog
import _mpl_imports as mpl
from yt.funcs import \
     ensure_list, \
     get_image_suffix

def get_canvas(name):
    suffix = get_image_suffix(name)
    
    if suffix == '':
        suffix = '.png'
    if suffix == ".png":
        canvas_cls = mpl.FigureCanvasAgg
    elif suffix == ".pdf":
        canvas_cls = mpl.FigureCanvasPdf
    elif suffix in (".eps", ".ps"):
        canvas_cls = mpl.FigureCanvasPS
    else:
        mylog.warning("Unknown suffix %s, defaulting to Agg", suffix)
        canvas_cls = FigureCanvasAgg
    return canvas_cls

def invalidate_plot(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        rv = f(*args, **kwargs)
        args[0]._plot_valid = False
        return rv
    return newfunc

class FigureContainer(dict):
    def __init__(self):
        super(dict, self).__init__()

    def __missing__(self, key):
        figure = mpl.matplotlib.figure.Figure((10, 8))
        self[key] = figure
        return self[key]

class AxesContainer(dict):
    def __init__(self, fig_container):
        self.fig_container = fig_container
        super(dict, self).__init__()

    def __missing__(self, key):
        figure = self.fig_container[key]
        self[key] = figure.add_subplot(111)
        return self[key]

class ProfilePlot(object):
    r"""
    Create a 1d profile plot from a data source or from a list 
    of profile objects.

    Given a data object (all_data, region, sphere, etc.), an x field, 
    and a y field (or fields), this will create a one-dimensional profile 
    of the average (or total) value of the y field in bins of the x field.

    This can be used to create profiles from given fields or to plot 
    multiple profiles created from 
    `yt.data_objects.profiles.create_profile`.
    
    Parameters
    ----------
    data_source : AMR4DData Object
        The data object to be profiled, such as all_data, region, or 
        sphere.
    x_field : str
        The binning field for the profile.
    y_fields : str or list
        The field or fields to be profiled.
    weight_field : str
        The weight field for calculating weighted averages.  If None, 
        the profile values are the sum of the field values within the bin.
        Otherwise, the values are a weighted average.
        Default : "CellMassMsun".
    n_bins : int
        The number of bins in the profile.
        Default: 64.
    accumulation : bool
        If True, the profile values for a bin N are the cumulative sum of 
        all the values from bin 0 to N.
        Default: False.
    fractional : If True the profile values are divided by the sum of all 
        the profile data such that the profile represents a probability 
        distribution function.
    label : str or list of strings
        If a string, the label to be put on the line plotted.  If a list, 
        this should be a list of labels for each profile to be overplotted.
        Default: None.
    plot_spec : dict or list of dicts
        A dictionary or list of dictionaries containing plot keyword 
        arguments.  For example, dict(color="red", linestyle=":").
        Default: None.
    profiles : list of profiles
        If not None, a list of profile objects created with 
        `yt.data_objects.profiles.create_profile`.
        Default: None.

    Examples
    --------

    This creates profiles of a single dataset.

    >>> pf = load("enzo_tiny_cosmology/DD0046/DD0046")
    >>> ad = pf.h.all_data()
    >>> plot = ProfilePlot(ad, "Density", ["Temperature", "x-velocity"], 
                           weight_field="CellMassMsun",
                           plot_spec=dict(color='red', linestyle="--"))
    >>> plot.save()

    This creates profiles from a time series object.
    
    >>> es = simulation("AMRCosmology.enzo", "Enzo")
    >>> es.get_time_series()

    >>> profiles = []
    >>> labels = []
    >>> plot_specs = []
    >>> for pf in es[-4:]:
    ...     ad = pf.h.all_data()
    ...     profiles.append(create_profile(ad, ["Density"],
    ...                                    fields=["Temperature",
    ...                                            "x-velocity"]))
    ...     labels.append(pf.current_redshift)
    ...     plot_specs.append(dict(linestyle="--", alpha=0.7))
    >>>
    >>> plot = ProfilePlot.from_profiles(profiles, labels=labels,
    ...                                  plot_specs=plot_specs)
    >>> plot.save()

    Use plot_line_property to change line properties of one or all profiles.
    
    """
    x_log = None
    y_log = None
    z_log = None
    x_title = None
    y_title = None

    _plot_valid = False

    def __init__(self, data_source, x_field, y_fields, 
                 weight_field="CellMassMsun", n_bins=64, 
                 accumulation=False, fractional=False,
                 label=None, plot_spec=None, profiles=None):
        self.y_log = {}
        self.y_title = {}
        if profiles is None:
            self.profiles = [create_profile(data_source, [x_field], n_bins,
                                            fields=ensure_list(y_fields),
                                            weight_field=weight_field,
                                            accumulation=accumulation,
                                            fractional=fractional)]
        else:
            self.profiles = ensure_list(profiles)
        
        self.label = label
        if not isinstance(self.label, list):
            self.label = [self.label] * len(self.profiles)

        self.plot_spec = plot_spec
        if self.plot_spec is None:
            self.plot_spec = [dict() for p in self.profiles]
        if not isinstance(self.plot_spec, list):
            self.plot_spec = [self.plot_spec.copy() for p in self.profiles]
        
        self._setup_plots()
        
    def save(self, name=None):
        r"""
        Saves a 1d profile plot.

        Parameters
        ----------
        name : str
            The output file keyword.
        
        """
        if not self._plot_valid: self._setup_plots()
        unique = set(self.figures.values())
        if len(unique) < len(self.figures):
            figiter = izip(xrange(len(unique)), sorted(unique))
        else:
            iters = self.figures.iteritems()
        if name is None:
            if len(self.profiles) == 1:
                prefix = self.profiles[0].pf
            else:
                prefix = "Multi-data"
            name = "%s.png" % prefix
        suffix = get_image_suffix(name)
        prefix = name[:name.rfind(suffix)]
        if not suffix:
            suffix = ".png"
        canvas_cls = get_canvas(name)
        for uid, fig in iters:
            canvas = canvas_cls(fig)
            fn = "%s_1d-Profile_%s_%s%s" % \
              (prefix, self.profiles[0].x_field, uid, suffix)
            mylog.info("Saving %s", fn)
            canvas.print_figure(fn)
        return self

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


    def _repr_html_(self):
        """Return an html representation of the plot object. Will display as a
        png for each WindowPlotMPL instance in self.plots"""
        ret = ''
        unique = set(self.figures.values())
        if len(unique) < len(self.figures):
            figiter = izip(xrange(len(unique)), sorted(unique))
        else:
            iters = self.figures.iteritems()
        for uid, fig in iters:
            canvas = mpl.FigureCanvasAgg(fig)
            f = cStringIO.StringIO()
            canvas.print_figure(f)
            f.seek(0)
            img = base64.b64encode(f.read())
            ret += '<img src="data:image/png;base64,%s"><br>' % img
        return ret

    def _setup_plots(self):
        self.figures = FigureContainer()
        self.axes = AxesContainer(self.figures)
        for i, profile in enumerate(self.profiles):
            for field, field_data in profile.field_data.items():
                self.axes[field].plot(profile.x[:-1], field_data, 
                                      label=self.label[i],
                                      **self.plot_spec[i])
        
        # This relies on 'profile' leaking
        for fname, axes in self.axes.items():
            xscale, yscale = self._get_field_log(fname, profile)
            xtitle, ytitle = self._get_field_title(fname, profile)
            axes.set_xscale(xscale)
            axes.set_yscale(yscale)
            axes.set_xlabel(xtitle)
            axes.set_ylabel(ytitle)
            axes.legend(loc="best")
        self._plot_valid = True

    @classmethod
    def from_profiles(cls, profiles, labels=None, plot_specs=None):
        r"""
        Instantiate a ProfilePlot object from a list of profiles 
        created with `yt.data_objects.profiles.create_profile`.

        Parameters
        ----------
        profiles : list of profiles
            If not None, a list of profile objects created with 
            `yt.data_objects.profiles.create_profile`.
        labels : list of strings
            A list of labels for each profile to be overplotted.
            Default: None.
        plot_specs : list of dicts
            A list of dictionaries containing plot keyword 
            arguments.  For example, [dict(color="red", linestyle=":")].
            Default: None.

        Examples
        --------

        >>> es = simulation("AMRCosmology.enzo", "Enzo")
        >>> es.get_time_series()

        >>> profiles = []
        >>> labels = []
        >>> plot_specs = []
        >>> for pf in es[-4:]:
        ...     ad = pf.h.all_data()
        ...     profiles.append(create_profile(ad, ["Density"],
        ...                                    fields=["Temperature",
        ...                                            "x-velocity"]))
        ...     labels.append(pf.current_redshift)
        ...     plot_specs.append(dict(linestyle="--", alpha=0.7))
        >>>
        >>> plot = ProfilePlot.from_profiles(profiles, labels=labels,
        ...                                  plot_specs=plot_specs)
        >>> plot.save()
        
        """
        if labels is not None and len(profiles) != len(labels):
            raise RuntimeError("Profiles list and labels list must be the same size.")
        if plot_specs is not None and len(plot_specs) != len(profiles):
            raise RuntimeError("Profiles list and plot_specs list must be the same size.")
        obj = cls(None, None, None, profiles=profiles, label=labels,
                  plot_spec=plot_specs)
        return obj

    @invalidate_plot
    def set_line_property(self, property, value, index=None):
        r"""
        Set properties for one or all lines to be plotted.

        Parameters
        ----------
        property : str
            The line property to be set.
        value : str, int, float
            The value to set for the line property.
        index : int
            The index of the profile in the list of profiles to be 
            changed.  If None, change all plotted lines.
            Default : None.

        Examples
        --------

        Change all the lines in a plot
        plot.set_line_property("linestyle", "-")

        Change a single line.
        plot.set_line_property("linewidth", 4, index=0)
        
        """
        if index is None:
            specs = self.plot_spec
        else:
            specs = [self.plot_spec[index]]
        for spec in specs:
            spec[property] = value
        return self
            
    def _get_field_log(self, field_y, profile):
        pf = profile.data_source.pf
        yfi = pf.field_info[field_y]
        if self.x_log is None:
            x_log = profile.x_log
        else:
            x_log = self.x_log
        if field_y in self.y_log:
            y_log = self.y_log[field_y]
        else:
            y_log = yfi.take_log
        scales = {True: 'log', False: 'linear'}
        return scales[x_log], scales[y_log]

    def _get_field_label(self, field, field_info):
        units = field_info.get_units()
        field_name = field_info.display_name
        if field_name is None:
            field_name = r'$\rm{'+field+r'}$'
        elif field_name.find('$') == -1:
            field_name = r'$\rm{'+field+r'}$'
        if units is None or units == '':
            label = field_name
        else:
            label = field_name+r'$\/\/('+units+r')$'
        return label

    def _get_field_title(self, field_y, profile):
        pf = profile.data_source.pf
        field_x = profile.x_field
        xfi = pf.field_info[field_x]
        yfi = pf.field_info[field_y]
        x_title = self.x_title or self._get_field_label(field_x, xfi)
        y_title = self.y_title.get(field_y, None) or \
                    self._get_field_label(field_y, yfi)
        return (x_title, y_title)
            

class PhasePlot(ImagePlotContainer):
    r"""
    Create a 2d profile (phase) plot from a data source or from 
    profile object created with 
    `yt.data_objects.profiles.create_profile`.

    Given a data object (all_data, region, sphere, etc.), an x field, 
    y field, and z field (or fields), this will create a two-dimensional 
    profile of the average (or total) value of the z field in bins of the 
    x and y fields.
    
    Parameters
    ----------
    data_source : AMR3DData Object
        The data object to be profiled, such as all_data, region, or 
        sphere.
    x_field : str
        The x binning field for the profile.
    y_field : str
        The y binning field for the profile.
    z_fields : str or list
        The field or fields to be profiled.
    weight_field : str
        The weight field for calculating weighted averages.  If None, 
        the profile values are the sum of the field values within the bin.
        Otherwise, the values are a weighted average.
        Default : "CellMassMsun".
    x_bins : int
        The number of bins in x field for the profile.
        Default: 128.
    y_bins : int
        The number of bins in y field for the profile.
        Default: 128.
    accumulation : bool or list of bools
        If True, the profile values for a bin n are the cumulative sum of 
        all the values from bin 0 to n.  If -True, the sum is reversed so 
        that the value for bin n is the cumulative sum from bin N (total bins) 
        to n.  A list of values can be given to control the summation in each
        dimension independently.
        Default: False.
    fractional : If True the profile values are divided by the sum of all 
        the profile data such that the profile represents a probability 
        distribution function.
    profile : profile object
        If not None, a profile object created with 
        `yt.data_objects.profiles.create_profile`.
        Default: None.
    fontsize: int
        Font size for all text in the plot.
        Default: 18.
    font_color : str
        Color for all text in the plot.
        Default: "black"
    figure_size : int
        Size in inches of the image.
        Default: 10 (10x10)

    Examples
    --------

    >>> pf = load("enzo_tiny_cosmology/DD0046/DD0046")
    >>> ad = pf.h.all_data()
    >>> plot = PhasePlot(ad, "Density", "Temperature", ["CellMassMsun"],
                         weight_field=None)
    >>> plot.save()

    >>> # Change plot properties.
    >>> plot.set_cmap("CellMassMsun", "jet")
    >>> plot.set_zlim("CellMassMsun", 1e8, 1e13)
    >>> plot.set_title("CellMassMsun", "This is a phase plot")
    
    """
    x_log = None
    y_log = None
    x_title = None
    y_title = None
    z_title = None
    plot_title = None
    _plot_valid = False
    _plot_type = 'Phase'

    def __init__(self, data_source, x_field, y_field, z_fields,
                 weight_field="CellMassMsun", x_bins=128, y_bins=128,
                 accumulation=False, fractional=False,
                 profile=None, fontsize=18, font_color="black", figure_size=8.0):
        self.plot_title = {}
        self.z_log = {}
        self.z_title = {}

        if profile is None:
            profile = create_profile(data_source,
               [x_field, y_field], [x_bins, y_bins],
               fields = ensure_list(z_fields),
               weight_field = weight_field,
               accumulation=accumulation,
               fractional=fractional)
        self.profile = profile
        ImagePlotContainer.__init__(self, data_source, profile.field_data.keys(),
                                    figure_size, fontsize)
        # This is a fallback, in case we forget.
        self._setup_plots()

    def _get_field_title(self, field_z, profile):
        pf = profile.data_source.pf
        field_x = profile.x_field
        field_y = profile.y_field
        xfi = pf.field_info[field_x]
        yfi = pf.field_info[field_y]
        zfi = pf.field_info[field_z]
        x_title = self.x_title or self._get_field_label(field_x, xfi)
        y_title = self.y_title or self._get_field_label(field_y, yfi)
        z_title = self.z_title.get(field_z, None) or \
                    self._get_field_label(field_z, zfi)
        return (x_title, y_title, z_title)

    def _get_field_label(self, field, field_info):
        units = field_info.get_units()
        field_name = field_info.display_name
        if field_name is None:
            field_name = r'$\rm{'+field+r'}$'
        elif field_name.find('$') == -1:
            field_name = r'$\rm{'+field+r'}$'
        if units is None or units == '':
            label = field_name
        else:
            label = field_name+r'$\/\/('+units+r')$'
        return label
        
    def _get_field_log(self, field_z, profile):
        pf = profile.data_source.pf
        zfi = pf.field_info[field_z]
        if self.x_log is None:
            x_log = profile.x_log
        else:
            x_log = self.x_log
        if self.y_log is None:
            y_log = profile.y_log
        else:
            y_log = self.y_log
        if field_z in self.z_log:
            z_log = self.z_log[field_z]
        else:
            z_log = zfi.take_log
        scales = {True: 'log', False: 'linear'}
        return scales[x_log], scales[y_log], scales[z_log]

    def _setup_plots(self):
        for f, data in self.profile.field_data.items():
            fig = None
            axes = None
            cax = None
            if f in self.plots:
                if self.plots[f].figure is not None:
                    fig = self.plots[f].figure
                    axes = self.plots[f].axes
                    cax = self.plots[f].cax

            size = (self.figure_size, self.figure_size)
            x_scale, y_scale, z_scale = self._get_field_log(f, self.profile)
            x_title, y_title, z_title = self._get_field_title(f, self.profile)
            if f in self.plots:
                zlim = [self.plots[f].zmin, self.plots[f].zmax]
            else:
                if z_scale == 'log':
                    zmin = data[data > 0.0].min()
                else:
                    zmin = data.min()
                zlim = [zmin, data.max()]                
            
            fp = self._font_properties
            self.plots[f] = PhasePlotMPL(self.profile.x, self.profile.y, data, 
                                         x_scale, y_scale, z_scale,
                                         self._colormaps[f], zlim, size, fp.get_size(),
                                         fig, axes, cax)
            self.plots[f].axes.xaxis.set_label_text(x_title)
            self.plots[f].axes.yaxis.set_label_text(y_title)
            self.plots[f].cax.yaxis.set_label_text(z_title)
            if f in self.plot_title:
                self.plots[f].axes.set_title(self.plot_title[f])

            if self._font_color is not None:
                ax = self.plots[f].axes
                cbax = self.plots[f].cb.ax
                labels = \
                  ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels() + \
                  cbax.yaxis.get_ticklabels() + \
                  [ax.xaxis.label, ax.yaxis.label, cbax.yaxis.label]
                for label in labels:
                    label.set_color(self._font_color)

    def save(self, name=None, mpl_kwargs=None):
        r"""
        Saves a 2d profile plot.

        Parameters
        ----------
        name : str
            The output file keyword.
        mpl_kwargs : dict
           A dict of keyword arguments to be passed to matplotlib.

        >>> plot.save(mpl_kwargs={'bbox_inches':'tight'})
        
        """

        if not self._plot_valid: self._setup_plots()
        if mpl_kwargs is None: mpl_kwargs = {}
        for f in self.profile.field_data:
            middle = "2d-Profile_%s_%s_%s" % (self.profile.x_field, 
                                              self.profile.y_field, f)
            if name is None:
                prefix = self.profile.pf
                name = "%s.png" % prefix
            suffix = get_image_suffix(name)
            prefix = name[:name.rfind(suffix)]
            fn = "%s_%s%s" % (prefix, middle, suffix)
            if not suffix:
                suffix = ".png"
            self.plots[f].save(fn, mpl_kwargs)

    @property
    def fields(self):
        return self.plots.keys()

    @invalidate_plot
    def set_title(self, field, title):
        """Set a title for the plot.

        Parameters
        ----------
        field : str
            The z field of the plot to add the title.
        title : str
            The title to add.

        Examples
        --------

        >>> plot.set_title("CellMassMsun", "This is a phase plot")
        
        """

        self.plot_title[field] = title

    @invalidate_plot
    def reset_plot(self):
        self.plots = {}
        return self
            
    def run_callbacks(self, *args):
        raise NotImplementedError
    def setup_callbacks(self, *args):
        raise NotImplementedError

class PhasePlotMPL(WindowPlotMPL):
    def __init__(self, x_data, y_data, data, 
                 x_scale, y_scale, z_scale, cmap,
                 zlim, size, fontsize, figure, axes, cax):
        self._draw_colorbar = True
        self._draw_axes = True
        self._cache_layout(size, fontsize)

        # Make room for a colorbar
        self.input_size = size
        self.fsize = [size[0] + self._cbar_inches[self._draw_colorbar], size[1]]

        # Compute layout
        axrect, caxrect = self._get_best_layout(fontsize=fontsize)
        if np.any(np.array(axrect) < 0):
            mylog.warning('The axis ratio of the requested plot is very narrow.  '
                          'There is a good chance the plot will not look very good, '
                          'consider making the plot manually using FixedResolutionBuffer '
                          'and matplotlib.')
            axrect  = (0.07, 0.10, 0.80, 0.80)
            caxrect = (0.87, 0.10, 0.04, 0.80)
        ImagePlotMPL.__init__(self, self.fsize, axrect, caxrect, zlim,
                              figure, axes, cax)
        self._init_image(x_data, y_data, data, x_scale, y_scale, z_scale,
                         zlim, cmap)

    def _init_image(self, x_data, y_data, image_data, 
                    x_scale, y_scale, z_scale, zlim, cmap):
        """Store output of imshow in image variable"""
        
        if (z_scale == 'log'):
            norm = matplotlib.colors.LogNorm(zlim[0], zlim[1])
        elif (z_scale == 'linear'):
            norm = matplotlib.colors.Normalize(zlim[0], zlim[1])
        self.image = None
        self.cb = None
        self.image = self.axes.pcolormesh(x_data, y_data, image_data,
                                          norm=norm, 
                                          cmap=cmap)
        self.axes.set_xscale(x_scale)
        self.axes.set_yscale(y_scale)
        self.cb = self.figure.colorbar(self.image, self.cax)
        if z_scale == 'linear':
            self.cb.formatter.set_scientific(True)
            self.cb.formatter.set_powerlimits((-2,3))
            self.cb.update_ticks()
