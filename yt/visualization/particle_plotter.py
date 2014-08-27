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
import types

from functools import wraps
from itertools import izip
import matplotlib
import numpy as np
import cStringIO

from yt.utilities.exceptions import \
    YTNotInsideNotebook
from yt.utilities.logger import ytLogger as mylog
import _mpl_imports as mpl
from yt.funcs import \
    ensure_list, \
    get_image_suffix, \
    get_ipython_api_version
from yt.units.unit_object import Unit
from .profile_plotter import \
    get_canvas, \
    invalidate_plot, \
    sanitize_label

class ParticlePlot(object):
    r"""
    Create a particle scatter plot from a data source.

    Given a data object (all_data, region, sphere, etc.), an x field, 
    and a y field (both of particle type), this will create a scatter
    plot with one marker for each particle.

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
        arguments.  This will be passed For example, dict('c'='r', 'marker'='.').
        Default: dict('c'='b', 'marker'='.', 'linestyle'='None', 'markersize'=8)

    Examples
    --------

    >>> import yt
    >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    >>> ad = ds.all_data()
    >>> plot = yt.ParticlePlot(ad, 'particle_position_x', 'particle_velocity_x')
    >>> plot.save()

    Use set_line_property to change line properties.
    
    """
    x_log = None
    y_log = None
    z_log = None
    x_title = None
    y_title = None
    x_lim = (None, None)
    y_lim = (None, None)
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
        if name is None:
            prefix = self.data_source.ds
            name = "%s.png" % prefix
        suffix = get_image_suffix(name)
        prefix = name[:name.rfind(suffix)]
        xfn = self.x_field
        if isinstance(xfn, types.TupleType):
            xfn = xfn[1]
        yfn = self.y_field
        if isinstance(yfn, types.TupleType):
            yfn = yfn[1]
        if not suffix:
            suffix = ".png"
        canvas_cls = get_canvas(name)
        canvas = canvas_cls(self.figure)
        fn = "%s_ScatterPlot_%s_%s%s" % (prefix, xfn, yfn, suffix)
        mylog.info("Saving %s", fn)
        canvas.print_figure(fn)
        return fn

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
        self.axis.cla()
        self.axis.plot(np.array(self.x_data), np.array(self.y_data),
                       **self.plot_spec)

        xscale, yscale = self._get_axis_log()
        xtitle, ytitle = self._get_axis_titles()

        self.axis.set_xscale(xscale)
        self.axis.set_yscale(yscale)

        self.axis.set_xlabel(xtitle)
        self.axis.set_ylabel(ytitle)

        self.axis.set_xlim(*self.x_lim)
        self.axis.set_ylim(*self.y_lim)

        self._plot_valid = True

    @invalidate_plot
    def set_line_property(self, property, value):
        r"""
        Set properties for one or all lines to be plotted.

        Parameters
        ----------
        property : str
            The line property to be set.
        value : str, int, float
            The value to set for the line property.

        Examples
        --------

        Change all the lines in a plot
        plot.set_line_property("linestyle", "-")

        Change a single line.
        plot.set_line_property("linewidth", 4, index=0)
        
        """
        specs = self.plot_spec
        specs[property] = value
        return self

    @invalidate_plot
    def set_xlog(self, log):
        """set the x-axis to log or linear.

        Parameters
        ----------
        field : string
            the field to set a transform
        log : boolean
            Log on/off.
        """
        self.x_log = log
        return self

    @invalidate_plot
    def set_ylog(self, log):
        """set the y-axis to log or linear.

        Parameters
        ----------
        field : string
            the field to set a transform
        log : boolean
            Log on/off.
        """
        self.y_log = log
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
        if field == self.x_field:
            self.x_data.convert_to_units(unit)
        elif field == self.y_field:
            self.y_data.convert_to_units(unit)
        else:
            raise KeyError("Field %s not in profile plot!" % (field))
        return self

    @invalidate_plot
    def set_xlim(self, xmin=None, xmax=None):
        """Sets the limits of the x field

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
        >>> pp = yt.ParticlePlot(ds.all_data(), 'particle_position_x', 'particle_velocity_x')
        >>> pp.set_xlim(0.1, 0.9)
        >>> pp.save()

        """
        self.x_lim = (xmin, xmax)
        return self

    @invalidate_plot
    def set_ylim(self, ymin=None, ymax=None):
        """Sets the limits for the y-axis of the plot.

        Parameters
        ----------

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
        >>> pp = yt.ParticlePlot(ds.all_data(), 'particle_position_x', 'particle_velocity_x')
        >>> pp.set_ylim(1e1, 1e8)
        >>> pp.save()

        """
        self.y_lim = (ymin, ymax)
        return self

    def _get_axis_log(self):

        xf, = self.data_source._determine_fields([self.x_field])
        xfi = self.data_source.ds._get_field_info(*xf)
        if self.x_log is None:
            x_log = xfi.take_log
        else:
            x_log = self.x_log

        yf, = self.data_source._determine_fields([self.y_field])
        yfi = self.data_source.ds._get_field_info(*yf)
        if self.y_log is None:
            y_log = yfi.take_log
        else:
            y_log = self.y_log
        
        scales = {True: 'log', False: 'linear'}
        return scales[x_log], scales[y_log]

    def _get_field_label(self, field, field_info, field_unit):
        field_unit = field_unit.latex_representation()
        field_name = field_info.display_name
        if isinstance(field, tuple): field = field[1]
        if field_name is None:
            field_name = r'$\rm{'+field+r'}$'
            field_name = r'$\rm{'+field.replace('_','\/').title()+r'}$'
        elif field_name.find('$') == -1:
            field_name = field_name.replace(' ','\/')
            field_name = r'$\rm{'+field_name+r'}$'
        if field_unit is None or field_unit == '' or field_unit == '1':
            label = field_name
        else:
            label = field_name+r'$\/\/('+field_unit+r')$'
        return label

    def _get_axis_titles(self):

        xfi = self.data_source.ds._get_field_info(self.x_field)
        x_unit = Unit(self.x_data.units, registry=self.data_source.ds.unit_registry)
        x_title = self._get_field_label(self.x_field, xfi, x_unit)

        yfi = self.data_source.ds._get_field_info(self.y_field)
        y_unit = Unit(self.y_data.units, registry=self.data_source.ds.unit_registry)
        y_title = self._get_field_label(self.y_field, yfi, y_unit)

        return (x_title, y_title)
