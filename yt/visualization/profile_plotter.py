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

from functools import wraps
from itertools import izip, repeat
import numpy as np
import cStringIO

from .image_writer import \
    write_image, apply_colormap
from yt.data_objects.profiles import \
     create_profile
from yt.utilities.lib import \
    write_png_to_string
from yt.data_objects.profiles import \
    BinnedProfile1D, \
    BinnedProfile2D
from .tick_locators import LogLocator, LinearLocator
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
        Default : "CellMass".
    n_bins : int
        The number of bins in the profile.
        Default: 64.
    accumulation : bool
        If True, the profile values for a bin N are the cumulative sum of 
        all the values from bin 0 to N.
        Default: False.
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

    >>> pf = load("DD0046/DD0046")
    >>> ad = pf.h.all_data()
    >>> plot = ProfilePlot(ad, "Density", ["Temperature", "x-velocity"], 
                           weight_field="CellMass",
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
    plot_spec = None
    x_log = None
    y_log = None
    x_title = None
    y_title = None

    _plot_valid = False

    def __init__(self, data_source, x_field, y_fields, 
                 weight_field="CellMass", n_bins=64, accumulation=False,
                 label=None, plot_spec=None, profiles=None):
        self.y_log = {}
        self.y_title = {}
        if profiles is None:
            self.profiles = [create_profile(data_source, [x_field], n_bins,
                                            fields=ensure_list(y_fields),
                                            weight_field=weight_field)]
        else:
            self.profiles = ensure_list(profiles)

        if accumulation:
            for profile in self.profiles:
                for field in profile.field_data:
                    profile.field_data[field] = \
                      profile.field_data[field].cumsum()
        
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
            

class PhasePlot(ProfilePlot):
    r"""
    Create a 2d profile (phase) plot from a data source or from a list of
    profile objects.

    Given a data object (all_data, region, sphere, etc.), an x field, 
    and a y field (or fields), this will create a one-dimensional profile 
    of the average (or total) value of the y field in bins of the x field.

    This can be used to create profiles from given fields or to plot 
    multiple profiles created from 
    `yt.data_objects.profiles.create_profile`.
    
    Parameters
    ----------
    data_source : AMR3DData Object
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
        Default : "CellMass".
    n_bins : int
        The number of bins in the profile.
        Default: 64.
    accumulation : bool
        If True, the profile values for a bin N are the cumulative sum of 
        all the values from bin 0 to N.
        Default: False.
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

    >>> pf = load("DD0046/DD0046")
    >>> ad = pf.h.all_data()
    >>> plot = ProfilePlot(ad, "Density", ["Temperature", "x-velocity"], 
                           weight_field="CellMass",
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
    plot_spec = None
    x_log = None
    y_log = None
    z_log = None
    x_title = None
    y_title = None
    z_title = None

    def __init__(self, data_source, x_field, y_field, z_fields,
                 weight_field="CellMassMsun", x_bins=128, y_bins=128,
                 plot_spec = None, profile = None):
        self.z_log = {}
        self.z_title = {}
        if profile is None:
            # We only 
            profile = create_profile(data_source,
               [x_field, y_field], [x_bins, y_bins],
               fields = ensure_list(z_fields),
               weight_field = weight_field)
        self.profile = profile
        # This is a fallback, in case we forget.
        self.plot_spec = plot_spec
        if self.plot_spec is None:
            self.plot_spec = {}
        self._setup_plots()

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
        xax = AxisSpec()
        xax.title = self.profile.x_bin_field
        xax.bounds = (self.profile._x_bins[0],
                      self.profile._x_bins[-1])
        xax.scale = {True: 'log', False: 'linear'}[self.profile._x_log]
        xax.calculate_ticks()

        yax = AxisSpec()
        yax.title = self.profile.y_bin_field
        yax.bounds = (self.profile._y_bins[0],
                      self.profile._y_bins[-1])
        yax.scale = {True: 'log', False: 'linear'}[self.profile._y_log]
        yax.calculate_ticks()

        cbar = ColorbarSpec()
        cbar.title = self._current_field
        if self.scale == 'log':
            nz = (self.profile[self._current_field] > 0)
            mi = self.profile[self._current_field][nz].min()
        else:
            mi = self.profile[self._current_field].min()
        ma = self.profile[self._current_field].max()
        cbar.bounds = (mi, ma)
        cbar.cmap = 'algae'
        cbar.scale = self.scale
        cbar.calculate_ticks()

        self.plot = ImagePlotContainer()
        self.plot.image = self.profile[self._current_field].transpose()
        self.plot.x_spec = xax
        self.plot.y_spec = yax
        self.plot.cbar = cbar

    def save(self, name=None):
        if not self._plot_valid: self._setup_plots()
        # We'll go with a mesh here, even if it's inappropriate
        im = self.image
        if self.cbar.scale == 'log':
            norm = mpl.matplotlib.colors.LogNorm()
        else:
            norm = mpl.matplotlib.colors.Normalize()
        pcm = axes.pcolormesh(x_bins, y_bins, self.image, norm=norm,
                              shading='flat', cmap = self.cbar.cmap,
                              rasterized=True)
        if self.x_spec.scale == 'log': axes.set_xscale("log")
        if self.y_spec.scale == 'log': axes.set_yscale("log")
        if self.x_spec.title is not None:
            axes.set_xlabel(self.x_spec.title)
        if self.y_spec.title is not None:
            axes.set_ylabel(self.y_spec.title)
        if isinstance(place, types.StringTypes):
            canvas = mpl.FigureCanvasAgg(figure)
            canvas.print_figure(place)
        return figure, axes

