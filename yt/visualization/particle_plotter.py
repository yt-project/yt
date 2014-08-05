"""
This is a simple mechanism for interfacing with Particle plots



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
from yt.units.unit_object import Unit

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
        super(FigureContainer, self).__init__()

    def __missing__(self, key):
        figure = mpl.matplotlib.figure.Figure((10, 8))
        self[key] = figure
        return self[key]

class AxesContainer(dict):
    def __init__(self, fig_container):
        self.fig_container = fig_container
        self.ylim = {}
        super(AxesContainer, self).__init__()

    def __missing__(self, key):
        figure = self.fig_container[key]
        self[key] = figure.add_subplot(111)
        return self[key]

    def __setitem__(self, key, value):
        super(AxesContainer, self).__setitem__(key, value)
        self.ylim[key] = (None, None)

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

class ParticlePlot(object):
    r"""
    Create a particle scatter plot from a data source.

    Given a data object (all_data, region, sphere, etc.), an x field, 
    and a y field (or fields), this will a scatter plot with one marker
    for each particle.

    Parameters
    ----------
    data_source : AMR3DData Object
        The data object to be profiled, such as all_data, region, or 
        sphere.
    x_field : str
        The field to plot on the x-axis.
    y_fields : str
        The field to plot on the y-axis.
    plot_spec : dict or list of dicts
        A dictionary or list of dictionaries containing plot keyword 
        arguments.  For example, dict(color="blue", linestyle=".").
        Default: None.

    Examples
    --------

    This creates profiles of a single dataset.

    >>> import yt
    >>> ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
    >>> ad = ds.all_data()
    >>> plot = ProfilePlot(ad, "density", ["temperature", "velocity_x"],
    ...                    weight_field="cell_mass",
    ...                    plot_spec=dict(color='red', linestyle="--"))
    >>> plot.save()

    This creates profiles from a time series object.

    >>> es = yt.simulation("AMRCosmology.enzo", "Enzo")
    >>> es.get_time_series()

    >>> profiles = []
    >>> labels = []
    >>> plot_specs = []
    >>> for ds in es[-4:]:
    ...     ad = ds.all_data()
    ...     profiles.append(create_profile(ad, ["density"],
    ...                                    fields=["temperature",
    ...                                            "velocity_x"]))
    ...     labels.append(ds.current_redshift)
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

    def __init__(self, data_source, x_field, y_field,
                 plot_spec=None):

        if plot_spec is None:
            plot_spec = {'c':'b', 'marker':'.', 'linestyle':'None', 'markersize':8}

        self.data_source = data_source
        self.x_field = x_field
        self.y_field = y_field
        self.plot_spec = plot_spec

        self.x_data = self.data_source[x_field]
        self.y_data = self.data_source[y_field]
        
        self.figure = mpl.matplotlib.figure.Figure((10, 8))
        self.axis = self.figure.add_subplot(111)
        self._setup_plots()

    def save(self, name=None):
        r"""
         Saves the scatter plot to disk.

         Parameters
         ----------
         name : str
             The output file keyword.

         """
        if not self._plot_valid:
            self._setup_plots()
        unique = set(self.figures.values())
        if len(unique) < len(self.figures):
            iters = izip(xrange(len(unique)), sorted(unique))
        else:
            iters = self.figures.iteritems()
        if name is None:
            if len(self.profiles) == 1:
                prefix = self.profiles[0].ds
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
        fns = []
        for uid, fig in iters:
            if isinstance(uid, types.TupleType):
                uid = uid[1]
            canvas = canvas_cls(fig)
            fns.append("%s_1d-Profile_%s_%s%s" % (prefix, xfn, uid, suffix))
            mylog.info("Saving %s", fns[-1])
            canvas.print_figure(fns[-1])
        return fns

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

        >>> import yt
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>> pp = ProfilePlot(ds.all_data(), 'density', 'temperature')
        >>> pp.show()

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
        canvas = mpl.FigureCanvasAgg(self.figure)
        f = cStringIO.StringIO()
        canvas.print_figure(f)
        f.seek(0)
        img = base64.b64encode(f.read())
        ret += '<img src="data:image/png;base64,%s"><br>' % img
        return ret

    def _setup_plots(self):
        self.axis.plot(np.array(self.x_data), np.array(self.y_data),
                       **self.plot_spec)

        xscale = self._get_field_log(self.x_field)
        yscale = self._get_field_log(self.y_field)

        xtitle = self._get_field_title(self.x_field)
        ytitle = self._get_field_title(self.y_field)

        self.axis.set_xscale(xscale)
        self.axis.set_yscale(yscale)

        self.axis.set_xlabel(xtitle)
        self.axis.set_ylabel(ytitle)

        self._plot_valid = True

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

        Parameters
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

    @invalidate_plot
    def set_xlim(self, xmin=None, xmax=None):
        """Sets the limits of the bin field

        Parameters
        ----------
        
        xmin : float or None
          The new x minimum.  Defaults to None, which leaves the xmin
          unchanged.

        xmax : float or None
          The new x maximum.  Defaults to None, which leaves the xmax
          unchanged.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>> pp = yt.ProfilePlot(ds.all_data(), 'density', 'temperature')
        >>> pp.set_xlim(1e-29, 1e-24)
        >>> pp.save()

        """
        for i, p in enumerate(self.profiles):
            if xmin is None:
                xmi = p.x_bins.min()
            else:
                xmi = xmin
            if xmax is None:
                xma = p.x_bins.max()
            else:
                xma = xmax
            extrema = {p.x_field: ((xmi, str(p.x.units)), (xma, str(p.x.units)))}
            units = {p.x_field: str(p.x.units)}
            for field in p.field_map.values():
                units[field] = str(p.field_data[field].units)
            self.profiles[i] = \
                create_profile(p.data_source, p.x_field,
                               n_bins=len(p.x_bins)-2,
                               fields=p.field_map.values(),
                               weight_field=p.weight_field,
                               accumulation=p.accumulation,
                               fractional=p.fractional,
                               extrema=extrema, units=units)
        return self

    @invalidate_plot
    def set_ylim(self, field, ymin=None, ymax=None):
        """Sets the plot limits for the specified field we are binning.

        Parameters
        ----------

        field : string or field tuple

        The field that we want to adjust the plot limits for.
        
        ymin : float or None
          The new y minimum.  Defaults to None, which leaves the ymin
          unchanged.

        ymax : float or None
          The new y maximum.  Defaults to None, which leaves the ymax
          unchanged.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>> pp = yt.ProfilePlot(ds.all_data(), 'density', ['temperature', 'x-velocity'])
        >>> pp.set_ylim('temperature', 1e4, 1e6)
        >>> pp.save()

        """
        for i, p in enumerate(self.profiles):
            if field is 'all':
                fields = self.axes.keys()
            else:
                fields = ensure_list(field)
            for profile in self.profiles:
                for field in profile.data_source._determine_fields(fields):
                    if field in profile.field_map:
                        field = profile.field_map[field]
                    self.axes.ylim[field] = (ymin, ymax)
                    # Continue on to the next profile.
                    break
        return self

    def _get_field_log(self, field):
        ds = self.data_source.ds
        f, = self.data_source._determine_fields([field])
        fi = ds._get_field_info(*f)
        do_log = fi.take_log
        scales = {True: 'log', False: 'linear'}
        return scales[do_log]

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
        elif field_unit is None or field_unit == '' or field_unit == '1':
            label = field_name
        else:
            field_unit = field_unit.latex_representation()
            label = field_name+r'$\/\/('+field_unit+r')$'
        return label

    def _get_field_title(self, field):
        ds = self.data_source.ds
        f, = self.data_source._determine_fields([field])
        fi = ds._get_field_info(*f)
        field_unit = Unit(fi.units, registry=self.data_source.ds.unit_registry)
        title = self._get_field_label(field, fi, field_unit)
        return title
