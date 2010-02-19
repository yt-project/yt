"""
This is an interface to MatPlotLib <http://matplotlib.sf.net> to plot
irregularly shaped grids, with the presumption that at any point we could have
data that is "hidden" in deeper levels of refinement.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.funcs import *
from FixedResolution import *
from matplotlib.figure import Figure

from matplotlib.backends.backend_agg \
    import FigureCanvasAgg as FigureCanvas
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import get_cmap
from matplotlib.ticker import LogFormatter, ScalarFormatter

cb_norms = {True: LogNorm, False: Normalize}
cb_ticks = {True: LogFormatter, False: ScalarFormatter}

class CallbackRegistryHandler(object):
    def __init__(self, plot):
        self.plot = plot
        self._callbacks = {}

    def __getitem__(self, item):
        if item not in self._callbacks:
            raise KeyError(item)
        cb = self._callbacks[item]

        @wraps(cb)
        def get_wrapper(*args, **kwargs):
            cbo = cb(*args, **kwargs)
            return self.plot.add_callback(cbo)

        return get_wrapper

    def __setitem__(self, item, val):
        self._callbacks[item] = val

    def __delitem__(self, item):
        del self._callbacks[item]

    def __iter__(self):
        for k in sorted(self._callbacks): yield k

    def keys(self):
        return self._callbacks.keys()

class AMRPlot(object):
    _invalid_plot = True
    def __init__(self, pf, figure, axes, fig_size = (10,8)):
        self.pf = pf

        if figure is None: figure = Figure(fig_size)
        # We accept the fig_size as the figure argument
        if isinstance(figure, types.TupleType):
            figure = Figure(figure)
        self._figure = figure

        if axes is None: axes = self._figure.add_subplot(1,1,1)
        self._axes = axes

        self._callbacks = []
        self._setup_callback_registry()

        self._invalidate_plot()

    def _setup_callback_registry(self):
        from Callbacks import callback_registry
        self.modify = CallbackRegistryHandler(self)
        for c in callback_registry.values():
            if not hasattr(c, '_type_name'): continue
            self.modify[c._type_name] = c

    def _validate_plot(self):
        # How do we handle the colorbars?
        # Maybe via a default callback?
        mylog.debug("Re-validating plot")
        if self._invalid_plot:
            self._redraw_plot()
            self._run_callbacks()
            self._invalid_plot = False

    def _run_callbacks(self):
        # We sublcass to avoid the problem of replacing self.image,
        # which is a collection
        self._axes.patches = []
        self._axes.texts = []
        for cb in self._callbacks:
            cb(self)

    def _invalidate_plot(self):
        mylog.debug("Invalidating plot")
        self._invalid_plot = True

    def _redraw_plot(self):
        pass

    def _generate_prefix(self):
        pass

    def save_image(self, prefix, format="png", override=True, force_save=False):
        self._validate_plot()
        if not override: prefix = self._generate_prefix(prefix)
        fn = ".".join([prefix, format])
        canvas = FigureCanvas(self._figure)
        if force_save:
            canvas.print_figure(fn)
        else:
            only_on_root(canvas.print_figure, fn)
        return fn

    def add_callback(self, func):
        """
        Add *func* as a callback to this plot.  *func* will be called with this
        plot as its first argument every time the plot is redrawn.  Returns the
        id of the callback (for use with :meth:`remove_callback`.)
        """
        self._callbacks.append(func)
        self._invalidate_plot()
        return len(self._callbacks)-1
    
class AMRColorBar(object):
    _field_name = None
    _field_label = None
    _limits = None
    def __init__(self, colorbar = None):
        self.label_args = {}
        self._colorbar = colorbar
        self.set_cmap(ytcfg.get("raven","colormap"))

    def set_cmap(self, cmap):
        if isinstance(cmap, types.StringTypes):
            cmap = get_cmap(cmap)
        self._cmap = cmap

    def set_limits(self, mi, ma):
        self._limits = (mi, ma)

    def autoscale(self):
        self._limits = None

    def set_label(self, label):
        self._field_label = label

    def __call__(self, plot):
        if plot._field_name != self._field_name:
            # Reset the field label
            self._field_name = plot._field_name
            self._field_label = plot._get_label()
            # If the field changed, we might need to change
            # our log state
            take_log = plot.pf.field_info[self._field_name].take_log
            self.norm = cb_norms[take_log]()
            if self._colorbar is not None:
                self._colorbar.set_norm(self.norm)
                self._colorbar.formatter = cb_ticks[take_log]()

        # We pre-empty any existing norm
        plot.image.set_cmap(self._cmap)
        plot.image.set_norm(self.norm)

        if self._limits is None:
            # We assume this is called *after* the buffer is generated
            mi = na.nanmin(plot.buff[self._field_name])
            ma = na.nanmax(plot.buff[self._field_name])
            limits = (mi, ma)
            mylog.debug("Setting limits to %s", limits)
        else:
            limits = self._limits

        self.norm.autoscale(limits)
        
        if self._colorbar is not None:
            self._colorbar.set_cmap(self._cmap)
            mylog.debug("Setting Label: %s", self._field_label)
            self._colorbar.set_label(self._field_label, **self.label_args)
            mylog.debug("Setting clim: %s", limits)
            self._colorbar.set_clim(*limits)
            self._colorbar.update_bruteforce(plot.image)

class AMRFieldPlot(AMRPlot):
    _field_name = None
    _field_label = None
    _antialias = True
    _period = (0.0, 0.0)
    _frb_args = None
    def __init__(self, data_source, field,
                 figure = None, axes = None,
                 setup_colorbar = True, periodic = False):
        self.data_source = data_source
        self._field_name = field
        if periodic: self._setup_periodicity()
        if setup_colorbar is True:
            fig_size = (10,8)
        else:
            fig_size = (8, 8)

        AMRPlot.__init__(self, data_source.pf, figure, axes, fig_size)

        # Now some modifications to the figure and axes
        self._figure.subplots_adjust(
            hspace=0.0, wspace=0.0,
            bottom=0.0, top=1.0, left=0.0, right=1.0)

        self.set_width(1, 'unitary')
        self._initialize_image(setup_colorbar)

    def _initialize_image(self, setup_colorbar):
        # We just need a simple temporary array
        l, b, width, height = self._axes.bbox.bounds
        temparray = na.ones((width, height))
        self.image = \
            self._axes.imshow(temparray, interpolation='nearest',
                              aspect=1.0, picker=True, origin='lower')
        self._axes.set_xticks(())
        self._axes.set_yticks(())
        self._axes.set_ylabel("")
        self._axes.set_xlabel("")
        if setup_colorbar is True:
            cb = self._figure.colorbar(self.image,
                extend='neither', shrink=0.95)
            self.add_callback(AMRColorBar(cb))
        elif isinstance(setup_colorbar, AMRColorBar):
            self.add_callback(setup_colorbar)
        else:
            self.add_callback(AMRColorBar(None))
        self._color_control = self._callbacks[-1]

    def set_clim(self, mi, ma):
        self._color_control.set_limits(mi, ma)
        # We can avoid all the other callbacks by only re-running this one
        self._color_control(self)

    def set_label(self, label):
        self._color_control.set_label(label)
        # We can avoid all the other callbacks by only re-running this one
        self._color_control(self)

    def set_cmap(self, cmap):
        self._color_control.set_cmap(cmap)
        # We can avoid all the other callbacks by only re-running this one
        self._color_control(self)

    def autoscale(self):
        self._color_control.autoscale()
        # We can avoid all the other callbacks by only re-running this one
        self._color_control(self)

    def set_width(self, width, unit):
        if isinstance(unit, types.StringTypes):
            unit = self.data_source.pf[str(unit)]
        # width is in code units
        self.width = width / unit
        xax, yax = lagos.x_dict[self.data_source.axis], \
                   lagos.y_dict[self.data_source.axis]
        self.xlim = (self.data_source.center[xax] - self.width*0.5,
                     self.data_source.center[xax] + self.width*0.5)
        self.ylim = (self.data_source.center[yax] - self.width*0.5,
                     self.data_source.center[yax] + self.width*0.5)
        self._invalidate_plot()

    def _setup_periodicity(self):
        axis = self.data_source.axis
        DLE = self.data.pf["DomainLeftEdge"]
        DRE = self.data.pf["DomainRightEdge"]
        DD = float(periodic)*(DRE - DLE)
        if axis < 3:
            xax = lagos.x_dict[axis]
            yax = lagos.y_dict[axis]
            self.xmin = DLE[xax] - DD[xax]
            self.xmax = DRE[xax] + DD[xax]
            self.ymin = DLE[yax] - DD[yax]
            self.ymax = DRE[yax] + DD[yax]
            self._period = (DD[xax], DD[yax])
        else:
            # Not quite sure how to deal with this, particularly
            # in the Orion case.  Cutting planes are tricky.
            self.xmin = self.ymin = 0.0
            self.xmax = self.ymax = 1.0

    def _redraw_plot(self):
        self._update_buffer()
        if self.image._A.size != self.buff.buff_size:
            self._axes.clear()
            aspect = (self.ylim[1]-self.ylim[0])/(self.xlim[1]-self.xlim[0])
            self.image = \
                self._axes.imshow(self.buff[self._field_name],
                    interpolation='nearest', aspect=1.0,
                    picker=True, origin='lower')
        else:
            self.image.set_data(self.buff[self._field_name])

    def _update_buffer(self):
        x0, x1 = self.xlim
        y0, y1 = self.ylim
        l, b, width, height = self._axes.bbox.bounds
        args = ((x0,x1,y0,y1), (width, height), self._antialias)
        if self._frb_args == args: return
        self._frb_args = args
        mylog.debug("Updating buffer: %s", self._frb_args)
        self.buff = FixedResolutionBuffer(self.data_source,
                        *self._frb_args)
        
    def _get_label(self):
        field_name = self._field_name
        proj = "proj" in self._type_name and \
               self.data_source._weight is None
        data_label = self.pf.field_info[field_name].get_label(proj)
        return data_label

class AMRSlicePlot(AMRFieldPlot):
    _type_name = "slice"

class AMRProjectionPlot(AMRFieldPlot):
    _type_name = "proj"

class AMRCuttingPlanePlot(AMRFieldPlot):
    pass

class AMRProfilePlot(AMRPlot):
    pass

class AMRPhasePlot(AMRPlot):
    pass

class AMRParticlePlot(AMRPlot):
    pass
