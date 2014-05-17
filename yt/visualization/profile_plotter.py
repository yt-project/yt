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


import __builtin__
import base64
import os
import types

from functools import wraps
from itertools import izip
import matplotlib
import numpy as np
import cStringIO
import __builtin__


from .base_plot_types import ImagePlotMPL
from .plot_container import \
    ImagePlotContainer, \
    log_transform, linear_transform
from yt.data_objects.profiles import \
     create_profile
from yt.utilities.exceptions import \
     YTNotInsideNotebook
from yt.utilities.logger import ytLogger as mylog
import _mpl_imports as mpl
from yt.funcs import \
     ensure_list, \
     get_image_suffix, \
     get_ipython_api_version

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
        canvas_cls = mpl.FigureCanvasAgg
    return canvas_cls

def invalidate_plot(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        rv = f(*args, **kwargs)
        args[0]._plot_valid = False
        args[0]._setup_plots()
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

def sanitize_label(label, nprofiles):
    label = ensure_list(label)
    
    if len(label) == 1:
        label = label * nprofiles
    
    if len(label) != nprofiles:
        raise RuntimeError("Number of labels must match number of profiles")

    for l in label:
        if l is not None and not isinstance(l, basestring):
            raise RuntimeError("All labels must be None or a string")

    return label

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
        Default : "cell_mass".
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
    >>> plot = ProfilePlot(ad, "density", ["temperature", "velocity_x"],
                           weight_field="cell_mass",
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
    ...     profiles.append(create_profile(ad, ["density"],
    ...                                    fields=["temperature",
    ...                                            "velocity_x"]))
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
                 weight_field="cell_mass", n_bins=64,
                 accumulation=False, fractional=False,
                 label=None, plot_spec=None, profiles=None):
        self.y_log = {}
        self.y_title = {}
        self.x_log = None
        if profiles is None:
            self.profiles = [create_profile(data_source, [x_field],
                                            n_bins=[n_bins],
                                            fields=ensure_list(y_fields),
                                            weight_field=weight_field,
                                            accumulation=accumulation,
                                            fractional=fractional)]
        else:
            self.profiles = ensure_list(profiles)
        
        self.label = sanitize_label(label, len(self.profiles))

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
        if not self._plot_valid:
            self._setup_plots()
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
        xfn = self.profiles[0].x_field
        if isinstance(xfn, types.TupleType):
            xfn = xfn[1]
        if not suffix:
            suffix = ".png"
        canvas_cls = get_canvas(name)
        for uid, fig in iters:
            if isinstance(uid, types.TupleType):
                uid = uid[1]
            canvas = canvas_cls(fig)
            fn = "%s_1d-Profile_%s_%s%s" % \
              (prefix, xfn, uid, suffix)
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
            for field, field_data in profile.items():
                self.axes[field].plot(np.array(profile.x),
                                      np.array(field_data),
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
            if any(self.label):
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
        if field == "all":
            self.x_log = log
            for field in self.profiles[0].field_data.keys():
                self.y_log[field] = log
        else:
            field, = self.profiles[0].data_source._determine_fields([field])
            if field == self.profiles[0].x_field:
                self.x_log = log
            elif field in self.profiles[0].field_data:
                self.y_log[field] = log
            else:
                raise KeyError("Field %s not in profile plot!" % (field))
        return self

    @invalidate_plot
    def set_unit(self, field, unit):
        """Sets a new unit for the requested field

        parameters
        ----------
        field : string
           The name of the field that is to be changed.

        new_unit : string or Unit object
           The name of the new unit.
        """
        for profile in self.profiles:
            if field == profile.x_field[1]:
                profile.set_x_unit(unit)
            elif field in self.profiles[0].field_map:
                profile.set_field_unit(field, unit)
            else:
                raise KeyError("Field %s not in profile plot!" % (field))
        return self

    def _get_field_log(self, field_y, profile):
        pf = profile.data_source.pf
        yf, = profile.data_source._determine_fields([field_y])
        yfi = pf._get_field_info(*yf)
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

    def _get_field_label(self, field, field_info, field_unit, fractional=False):
        field_unit = field_unit.latex_representation()
        field_name = field_info.display_name
        if isinstance(field, tuple): field = field[1]
        if field_name is None:
            field_name = r'$\rm{'+field+r'}$'
            field_name = r'$\rm{'+field.replace('_','\/').title()+r'}$'
        elif field_name.find('$') == -1:
            field_name = field_name.replace(' ','\/')
            field_name = r'$\rm{'+field_name+r'}$'
        if fractional:
            label = field_name + r'$\rm{\/Probability\/Density}$'
        elif field_unit is None or field_unit == '':
            label = field_name
        else:
            label = field_name+r'$\/\/('+field_unit+r')$'
        return label

    def _get_field_title(self, field_y, profile):
        pf = profile.data_source.pf
        field_x = profile.x_field
        xf, yf = profile.data_source._determine_fields(
            [field_x, field_y])
        xfi = pf._get_field_info(*xf)
        yfi = pf._get_field_info(*yf)
        x_unit = profile.x.units
        y_unit = profile.field_units[field_y]
        fractional = profile.fractional
        x_title = self.x_title or self._get_field_label(field_x, xfi, x_unit)
        y_title = self.y_title.get(field_y, None) or \
                    self._get_field_label(field_y, yfi, y_unit, fractional)

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
        Default : "cell_mass".
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
        Default: 8 (8x8)

    Examples
    --------

    >>> pf = load("enzo_tiny_cosmology/DD0046/DD0046")
    >>> ad = pf.h.all_data()
    >>> plot = PhasePlot(ad, "density", "temperature", ["cell_mass"],
                         weight_field=None)
    >>> plot.save()

    >>> # Change plot properties.
    >>> plot.set_cmap("cell_mass", "jet")
    >>> plot.set_zlim("cell_mass", 1e8, 1e13)
    >>> plot.set_title("cell_mass", "This is a phase plot")
    
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
                 weight_field="cell_mass", x_bins=128, y_bins=128,
                 accumulation=False, fractional=False,
                 profile=None, fontsize=18, font_color="black", figure_size=8.0):
        self.plot_title = {}
        self.z_log = {}
        self.z_title = {}
        self._initfinished = False
        self.x_log = None
        self.y_log = None

        if profile is None:
            profile = create_profile(data_source,
               [x_field, y_field],
               ensure_list(z_fields),
               n_bins = [x_bins, y_bins],
               weight_field = weight_field,
               accumulation=accumulation,
               fractional=fractional)
        self.profile = profile
        super(PhasePlot, self).__init__(data_source, figure_size, fontsize)
        # This is a fallback, in case we forget.
        self._setup_plots()
        self._initfinished = True

    def _get_field_title(self, field_z, profile):
        pf = profile.data_source.pf
        field_x = profile.x_field
        field_y = profile.y_field
        xf, yf, zf = profile.data_source._determine_fields(
            [field_x, field_y, field_z])
        xfi = pf._get_field_info(*xf)
        yfi = pf._get_field_info(*yf)
        zfi = pf._get_field_info(*zf)
        x_unit = profile.x.units
        y_unit = profile.y.units
        z_unit = profile.field_units[field_z]
        fractional = profile.fractional
        x_title = self.x_title or self._get_field_label(field_x, xfi, x_unit)
        y_title = self.y_title or self._get_field_label(field_y, yfi, y_unit)
        z_title = self.z_title.get(field_z, None) or \
                    self._get_field_label(field_z, zfi, z_unit, fractional)
        return (x_title, y_title, z_title)

    def _get_field_label(self, field, field_info, field_unit, fractional=False):
        field_unit = field_unit.latex_representation()
        field_name = field_info.display_name
        if isinstance(field, tuple): field = field[1]
        if field_name is None:
            field_name = r'$\rm{'+field+r'}$'
            field_name = r'$\rm{'+field.replace('_','\/').title()+r'}$'
        elif field_name.find('$') == -1:
            field_name = field_name.replace(' ','\/')
            field_name = r'$\rm{'+field_name+r'}$'
        if fractional:
            label = field_name + r'$\rm{\/Probability\/Density}$'
        elif field_unit is None or field_unit is '':
            label = field_name
        else:
            label = field_name+r'$\/\/('+field_unit+r')$'
        return label
        
    def _get_field_log(self, field_z, profile):
        pf = profile.data_source.pf
        zf, = profile.data_source._determine_fields([field_z])
        zfi = pf._get_field_info(*zf)
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
        for f, data in self.profile.items():
            fig = None
            axes = None
            cax = None
            draw_colorbar = True
            draw_axes = True
            zlim = (None, None)
            if f in self.plots:
                draw_colorbar = self.plots[f]._draw_colorbar
                draw_axes = self.plots[f]._draw_axes
                zlim = (self.plots[f].zmin, self.plots[f].zmax)
                if self.plots[f].figure is not None:
                    fig = self.plots[f].figure
                    axes = self.plots[f].axes
                    cax = self.plots[f].cax

            x_scale, y_scale, z_scale = self._get_field_log(f, self.profile)
            x_title, y_title, z_title = self._get_field_title(f, self.profile)

            if zlim == (None, None):
                if z_scale == 'log':
                    zmin = data[data > 0.0].min()
                    self._field_transform[f] = log_transform
                else:
                    zmin = data.min()
                    self._field_transform[f] = linear_transform
                zlim = [zmin, data.max()]

            fp = self._font_properties
            f = self.profile.data_source._determine_fields(f)[0]

            self.plots[f] = PhasePlotMPL(self.profile.x, self.profile.y, data,
                                         x_scale, y_scale, z_scale,
                                         self._colormaps[f], np.array(zlim),
                                         self.figure_size, fp.get_size(),
                                         fig, axes, cax)

            self.plots[f]._toggle_axes(draw_axes)
            self.plots[f]._toggle_colorbar(draw_colorbar)

            self.plots[f].axes.xaxis.set_label_text(x_title, fontproperties=fp)
            self.plots[f].axes.yaxis.set_label_text(y_title, fontproperties=fp)
            self.plots[f].cax.yaxis.set_label_text(z_title, fontproperties=fp)

            if f in self.plot_title:
                self.plots[f].axes.set_title(self.plot_title[f])

            ax = self.plots[f].axes
            cbax = self.plots[f].cb.ax
            labels = ((ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels() +
                       cbax.yaxis.get_ticklabels()) +
                      [ax.xaxis.label, ax.yaxis.label, cbax.yaxis.label])
            for label in labels:
                label.set_fontproperties(fp)
                if self._font_color is not None:
                    label.set_color(self._font_color)
        self._plot_valid = True

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
        names = []
        if not self._plot_valid:
            self._setup_plots()
        if mpl_kwargs is None:
            mpl_kwargs = {}
        if name is None:
            name = str(self.profile.pf)
        xfn = self.profile.x_field
        yfn = self.profile.y_field
        if isinstance(xfn, types.TupleType):
            xfn = xfn[1]
        if isinstance(yfn, types.TupleType):
            yfn = yfn[1]
        for f in self.profile.field_data:
            _f = f
            if isinstance(f, types.TupleType):
                _f = _f[1]
            middle = "2d-Profile_%s_%s_%s" % (xfn, yfn, _f)
            splitname = os.path.split(name)
            if splitname[0] != '' and not os.path.isdir(splitname[0]):
                os.makedirs(splitname[0])
            if os.path.isdir(name) and name != str(self.profile.pf):
                prefix = name + (os.sep if name[-1] != os.sep else '')
                prefix += str(self.profile.pf)
            else:
                prefix = name
            suffix = get_image_suffix(name)
            if suffix != '':
                for k, v in self.plots.iteritems():
                    names.append(v.save(name, mpl_kwargs))
                return names
            fn = "%s_%s%s" % (prefix, middle, suffix)
            names.append(fn)
            self.plots[f].save(fn, mpl_kwargs)
        return names

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

        >>> plot.set_title("cell_mass", "This is a phase plot")
        
        """

        self.plot_title[field] = title

    @invalidate_plot
    def reset_plot(self):
        self.plots = {}
        return self

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
        if field == "all":
            self.x_log = log
            self.y_log = log
            for field in self.profile.field_data:
                self.z_log[field] = log
        else:
            if field == self.profile.x_field[1]:
                self.x_log = log
            elif field == self.profile.y_field[1]:
                self.y_log = log
            elif field in self.profile.field_map:
                self.z_log[self.profile.field_map[field]] = log
            else:
                raise KeyError("Field %s not in phase plot!" % (field))
        return self

    @invalidate_plot
    def set_unit(self, field, unit):
        """Sets a new unit for the requested field

        parameters
        ----------
        field : string
           The name of the field that is to be changed.

        new_unit : string or Unit object
           The name of the new unit.
        """
        fields = [fd[1] for fd in self.profile.field_data]
        if field == self.profile.x_field[1]:
            self.profile.set_x_unit(unit)
        elif field == self.profile.y_field[1]:
            self.profile.set_y_unit(unit)
        elif field in fields:
            self.profile.set_field_unit(field, unit)
        else:
            raise KeyError("Field %s not in phase plot!" % (field))
        return self

    def run_callbacks(self, *args):
        raise NotImplementedError
    def setup_callbacks(self, *args):
        raise NotImplementedError


class PhasePlotMPL(ImagePlotMPL):
    def __init__(self, x_data, y_data, data,
                 x_scale, y_scale, z_scale, cmap,
                 zlim, figure_size, fontsize, figure, axes, cax):
        self._initfinished = False
        self._draw_colorbar = True
        self._draw_axes = True
        self._figure_size = figure_size

        # Compute layout
        fontscale = float(fontsize) / 18.0
        if fontscale < 1.0:
            fontscale = np.sqrt(fontscale)

        self._cb_size = 0.0375*figure_size
        self._ax_text_size = [1.1*fontscale, 0.9*fontscale]
        self._top_buff_size = 0.30*fontscale
        self._aspect = 1.0

        size, axrect, caxrect = self._get_best_layout()

        super(PhasePlotMPL, self).__init__(size, axrect, caxrect, zlim,
                                           figure, axes, cax)

        self._init_image(x_data, y_data, data, x_scale, y_scale, z_scale,
                         zlim, cmap)

        self._initfinished = True

    def _init_image(self, x_data, y_data, image_data,
                    x_scale, y_scale, z_scale, zlim, cmap):
        """Store output of imshow in image variable"""
        if (z_scale == 'log'):
            norm = matplotlib.colors.LogNorm(zlim[0], zlim[1])
        elif (z_scale == 'linear'):
            norm = matplotlib.colors.Normalize(zlim[0], zlim[1])
        self.image = None
        self.cb = None
        self.image = self.axes.pcolormesh(np.array(x_data), 
                                          np.array(y_data), 
                                          np.array(image_data.T),
                                          norm=norm, 
                                          cmap=cmap)
        self.axes.set_xscale(x_scale)
        self.axes.set_yscale(y_scale)
        self.cb = self.figure.colorbar(self.image, self.cax)
        if z_scale == 'linear':
            self.cb.formatter.set_scientific(True)
            self.cb.formatter.set_powerlimits((-2,3))
            self.cb.update_ticks()
