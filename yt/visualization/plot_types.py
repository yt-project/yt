"""
This is an interface to MatPlotLib <http://matplotlib.sf.net> to plot
irregularly shaped grids, with the presumption that at any point we could have
data that is "hidden" in deeper levels of refinement.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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

import numpy as na

from yt.funcs import *
from _mpl_imports import *
import _MPL
from .plot_modifications import callback_registry
from yt.utilities.definitions import \
    x_dict, \
    y_dict, \
    axis_names
from .color_maps import yt_colormaps, is_colormap
from yt.utilities.exceptions import YTNoDataInObjectError

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

class RavenPlot(object):

    datalabel = None
    colorbar = None
    xlim = None
    ylim = None
    def __init__(self, data, fields, figure = None, axes=None, size=(10,8)):
        self.data = data
        self.fields = fields
        self.size = size
        self.set_autoscale(True)
        self.im = defaultdict(lambda: "")
        self["ParameterFile"] = "%s" % self.data.pf
        self.pf = self.data.pf
        self.axis_names = {}
        if not figure:
            self._figure = matplotlib.figure.Figure(size)
        else:
            self._figure = figure
        if not axes:
            self._axes = self._figure.add_subplot(1,1,1)
        else:
            self._axes = axes
        self._callbacks = []
        self._setup_callback_registry()

    def set_autoscale(self, val):
        self.do_autoscale = val

    def __getitem__(self, item):
        return self.data[item] # Should be returned in CGS

    def save_image(self, prefix, format="png", override=True, force_save=False,
                   figure_canvas = None):
        """
        Save this plot image.  Will generate a filename based on the *prefix*,
        *format*.  *override* will force no filename generation beyond the
        prefix.
        """
        self._redraw_image()
        if not override:
            self._generate_prefix(prefix)
            my_prefix = self.prefix
        else:
            my_prefix = prefix
        fn = ".".join([my_prefix, format])
        if figure_canvas is None:
            figure_canvas = FigureCanvasAgg
        canvas = figure_canvas(self._figure)
        if force_save:
            canvas.print_figure(fn)
        else:
            only_on_root(canvas.print_figure, fn)
        self["Type"] = self._type_name
        self["GeneratedAt"] = self.data.pf.unique_identifier
        return fn

    def open_in_webbrowser(self):
        from mplh5canvas.backend_h5canvas import FigureCanvasH5Canvas
        canv = FigureCanvasH5Canvas(self._figure)
        canv.draw()

    def save_to_pdf(self, f):
        self._redraw_image()
        canvas = FigureCanvasPdf(self._figure)
        original_figure_alpha = self._figure.patch.get_alpha()
        self._figure.patch.set_alpha(0.0)
        original_axes_alpha = []
        for ax in self._figure.axes:
            patch = ax.patch
            original_axes_alpha.append(patch.get_alpha())
            patch.set_alpha(0.0)

        canvas.print_pdf(f)

        self._figure.set_alpha(original_figure_alpha)
        for ax, alpha in zip(self._figure.axes,original_axes_alpha):
            ax.patch.set_alpha(alpha)

    def _redraw_image(self):
        pass

    def _generate_prefix(self):
        pass

    def set_xlim(self, xmin, xmax):
        """
        Set the x boundaries of this plot.
        """
        self.xlim = (xmin,xmax)

    def set_ylim(self, ymin, ymax):
        """
        Set the y boundaries of this plot.
        """
        self.ylim = (ymin,ymax)


    def set_zlim(self, zmin, zmax, dex=None, nticks=None, ticks=None, minmaxtick=False):
        """
        Set the z boundaries of this plot.

        Only ONE of the following options can be specified. If all 3 are
        specified, they will be used in the following precedence order:

        * ``ticks`` - a list of floating point numbers at which to put ticks
        * ``minmaxtick`` - display DEFAULT ticks with min & max also displayed
        * ``nticks`` - if ticks not specified, can automatically determine a
          number of ticks to be evenly spaced in log space
        """
        # This next call fixes some things, but is slower...
        self._redraw_image()
        self.set_autoscale(False)
        if (zmin in (None,'min')) or (zmax in (None,'max')):    
            imbuff = self._axes.images[-1]._A
            if zmin == 'min':
                zmin = na.nanmin(imbuff[na.nonzero(imbuff)])
                if dex is not None:
                    zmax = min(zmin*10**(dex),na.nanmax(imbuff))
            if zmax == 'max':
                zmax = na.nanmax(imbuff)
                if dex is not None:
                    zmin = max(zmax/(10**(dex)),na.nanmin(imbuff))
        if self.colorbar is not None:
            if ticks is not None:
                ticks = na.sort(ticks)
                self.colorbar.locator = matplotlib.ticker.FixedLocator(ticks)
                self.colorbar.formatter = matplotlib.ticker.FixedFormatter(["%0.2e" % (x) for x in ticks])
            elif minmaxtick:
                if self.log_field: 
                    ticks = na.array(self.colorbar._ticker()[1],dtype='float')
                    ticks = [zmin] + ticks.tolist() + [zmax]
                    self.colorbar.locator = matplotlib.ticker.FixedLocator(ticks)
                    self.colorbar.formatter = matplotlib.ticker.FixedFormatter(["%0.2e" % (x) for x in ticks])
                else:
                    mylog.error('Sorry, we do not support minmaxtick for linear fields.  It likely comes close by default')
            elif nticks is not None:
                if self.log_field:
                    lin = na.linspace(na.log10(zmin),na.log10(zmax),nticks)
                    self.colorbar.locator = matplotlib.ticker.FixedLocator(10**lin)
                    self.colorbar.formatter = matplotlib.ticker.FixedFormatter(["%0.2e" % (10**x) for x in lin])
                else: 
                    lin = na.linspace(zmin,zmax,nticks)
                    self.colorbar.locator = matplotlib.ticker.FixedLocator(lin)
                    self.colorbar.formatter = matplotlib.ticker.FixedFormatter(["%0.2e" % x for x in lin])

            else:
                if hasattr(self,'_old_locator'):
                    self.colorbar.locator = self._old_locator
                if hasattr(self,'_old_formatter'):
                    self.colorbar.formatter = self._old_formatter
        self.norm.autoscale(na.array([zmin,zmax], dtype='float64'))
        self.image.changed()
        if self.colorbar is not None:
            mpl_notify(self.image, self.colorbar)

    def set_cmap(self, cmap):
        """
        Change the colormap of this plot to *cmap*.
        """
        if isinstance(cmap, types.StringTypes):
            if str(cmap) in yt_colormaps:
                cmap = yt_colormaps[str(cmap)]
            elif hasattr(matplotlib.cm, cmap):
                cmap = getattr(matplotlib.cm, cmap)
        if not is_colormap(cmap) and cmap is not None:
            raise RuntimeError("Colormap '%s' does not exist!" % str(cmap))
        else:
            self.cmap = cmap
        
    def __setitem__(self, item, val):
        self.im[item] = val

    def add_callback(self, func):
        """
        Add *func* as a callback to this plot.  *func* will be called with this
        plot as its first argument every time the plot is redrawn.  Returns the
        id of the callback (for use with :meth:`remove_callback`.)
        """
        self._callbacks.append(func)
        return len(self._callbacks)-1

    def remove_callback(self, id):
        """
        Given an *id*, remove that index in the callbacks list.
        """
        self._callbacks[id] = lambda a: None

    def _run_callbacks(self):
        self._axes.patches = []
        self._axes.collections = []
        self._axes.texts = []
        for cb in self._callbacks:
            cb(self)

    def set_label(self, label):
        """
        Set the datalabel to *label*.  (This has different meanings based on
        the plot.)
        """
        self.datalabel = label
        if self.colorbar != None: self.colorbar.set_label(str(label))

    def setup_domain_edges(self, axis, periodic=False):
        DLE = self.data.pf.domain_left_edge
        DRE = self.data.pf.domain_right_edge
        DD = float(periodic)*(DRE - DLE)
        if axis < 3:
            xax = x_dict[axis]
            yax = y_dict[axis]
            self.xmin = DLE[xax] - DD[xax]
            self.xmax = DRE[xax] + DD[xax]
            self.ymin = DLE[yax] - DD[yax]
            self.ymax = DRE[yax] + DD[yax]
            self._period = (DD[xax], DD[yax])
            self._edges = ( (DLE[xax], DRE[xax]), (DLE[yax], DRE[yax]) )
        else:
            # Not quite sure how to deal with this, particularly
            # in the Orion case.  Cutting planes are tricky.
            self.xmin = self.ymin = 0.0
            self.xmax = self.ymax = 1.0

    def _setup_callback_registry(self):
        self.modify = CallbackRegistryHandler(self)
        for c in callback_registry.values():
            if not hasattr(c, '_type_name'): continue
            self.modify[c._type_name] = c

    def _pretty_name(self):
        width = self.im.get("Width", "NA")
        unit = self.im.get("Unit", "NA")
        field = self.axis_names.get("Z", self.axis_names.get("Field1"))
        if hasattr(self.data, "_data_source"):
            data = self.data._data_source
        else:
            data = self.data
        return "%s: %s (%s %s) %s" % (self._type_name,
            field, width, unit, data)

class VMPlot(RavenPlot):
    _antialias = True
    _period = (0.0, 0.0)
    _edges = True
    def __init__(self, data, field, figure = None, axes = None,
                 use_colorbar = True, size=None, periodic = False):
        fields = ['X', 'Y', field, 'X width', 'Y width']
        if not size:
            size = (10,8)
            if not use_colorbar: size=(8,8)
        RavenPlot.__init__(self, data, fields, figure, axes, size=size)
        self._figure.subplots_adjust(hspace=0, wspace=0, bottom=0.0,
                                    top=1.0, left=0.0, right=1.0)
        self.setup_domain_edges(self.data.axis, periodic)
        self.cmap = None
        self.label_kws = {}
        self.__setup_from_field(field)
        self.__init_temp_image(use_colorbar)

    def __setup_from_field(self, field):
        self.set_log_field(self.pf.field_info[field].take_log)
        self.axis_names["Z"] = field

    def set_log_field(self, val):
        if val:
            self.log_field = True
            self.norm = matplotlib.colors.LogNorm()
            ttype = matplotlib.ticker.LogFormatter
        else:
            self.log_field = False
            self.norm = matplotlib.colors.Normalize()
            ttype = matplotlib.ticker.ScalarFormatter
        if self.colorbar:
            self.colorbar.set_norm(self.norm)
            self.colorbar.formatter = ttype()

    def __init_temp_image(self, setup_colorbar):
        temparray = na.ones(self.size)
        self.image = \
            self._axes.imshow(temparray, interpolation='nearest',
                             norm = self.norm, aspect=1.0, picker=True,
                             origin='lower')
        self._axes.set_xticks(())
        self._axes.set_yticks(())
        self._axes.set_ylabel("")
        self._axes.set_xlabel("")
        if setup_colorbar:
            self.colorbar = self._figure.colorbar(self._axes.images[-1], \
                                                  extend='neither', \
                                                  shrink=0.95)
            self._old_locator = self.colorbar.locator
            self._old_formatter = self.colorbar.formatter
        else:
            self.colorbar = None
        self.set_width(1,'unitary')

    def _get_buff(self, width=None):
        x0, x1 = self.xlim
        y0, y1 = self.ylim
        if width is None:
            l, b, width, height = mpl_get_bounds(self._axes.bbox)
        else:
            height = width
        self.pix = (width,height)
        # 'px' == pixel x, or x in the plane of the slice
        # 'x' == actual x
        aa = int(self._antialias)
        if self._edges is True or \
           x0 < self._edges[0][0] or x1 > self._edges[0][1] or \
           y0 < self._edges[1][0] or y1 > self._edges[1][1]:
            check_period = 1
        else:
            check_period = 0
        buff = _MPL.Pixelize(self.data['px'],
                            self.data['py'],
                            self.data['pdx'],
                            self.data['pdy'],
                            self[self.axis_names["Z"]],
                            int(height), int(width),
                            (x0, x1, y0, y1),aa,self._period,
                            check_period).transpose()
        return buff

    def _redraw_image(self, *args):
        buff = self._get_buff()
        if self[self.axis_names["Z"]].size == 0:
            raise YTNoDataInObjectError(self.data)
        mylog.debug("Received buffer of min %s and max %s (data: %s %s)",
                    na.nanmin(buff), na.nanmax(buff),
                    self[self.axis_names["Z"]].min(),
                    self[self.axis_names["Z"]].max())
        if self.log_field:
            bI = na.where(buff > 0)
            if len(bI[0]) == 0:
                newmin = 1e-99
                newmax = 1e-99
            else:
                newmin = na.nanmin(buff[bI])
                newmax = na.nanmax(buff[bI])
        else:
            newmin = na.nanmin(buff)
            newmax = na.nanmax(buff)
        aspect = (self.ylim[1]-self.ylim[0])/(self.xlim[1]-self.xlim[0])
        if self.image._A.size != buff.size:
            self._axes.clear()
            self.image = \
                self._axes.imshow(buff, interpolation='nearest', norm = self.norm,
                                aspect=aspect, picker=True, origin='lower')
        else:
            self.image.set_data(buff)
        if self._axes.get_aspect() != aspect: self._axes.set_aspect(aspect)
        if self.do_autoscale:
            self.norm.autoscale(na.array((newmin,newmax), dtype='float64'))
        self._reset_image_parameters()
        self._run_callbacks()

    def _reset_image_parameters(self):
        self._axes.set_xticks(())
        self._axes.set_yticks(())
        self._axes.set_ylabel("")
        self._axes.set_xlabel("")
        if self.cmap:
            self.image.set_cmap(self.cmap)
        if self.colorbar != None:
            self.image.set_norm(self.norm)
            self.colorbar.set_norm(self.norm)
            if self.cmap: self.colorbar.set_cmap(self.cmap)
            if self.do_autoscale: mpl_notify(self.image, self.colorbar)
        self._autoset_label()

    def set_xlim(self, xmin, xmax):
        self.xlim = (xmin,xmax)

    def set_ylim(self, ymin, ymax):
        self.ylim = (ymin,ymax)

    def _generate_prefix(self, prefix):
        self.prefix = "_".join([prefix, self._type_name, \
            axis_names[self.data.axis], self.axis_names['Z']])
        self["Field1"] = self.axis_names["Z"]
        self["Field2"] = None
        self["Field3"] = None

    def set_width(self, width, unit):
        self["Unit"] = str(unit)
        self["Width"] = float(width)
        if isinstance(unit, types.StringTypes):
            unit = self.data.pf[str(unit)]
        self.width = width / unit
        self._refresh_display_width()

    def _refresh_display_width(self, width=None):
        if width:
            self.width = width
        else:
            width = self.width
        if iterable(width):
            width_x, width_y = width
        else:
            width_x = width
            width_y = width
        l_edge_x = self.data.center[x_dict[self.data.axis]] - width_x/2.0
        r_edge_x = self.data.center[x_dict[self.data.axis]] + width_x/2.0
        l_edge_y = self.data.center[y_dict[self.data.axis]] - width_y/2.0
        r_edge_y = self.data.center[y_dict[self.data.axis]] + width_y/2.0
        self.set_xlim(max(l_edge_x,self.xmin), min(r_edge_x,self.xmax))
        self.set_ylim(max(l_edge_y,self.ymin), min(r_edge_y,self.ymax))
        self._redraw_image()

    def autoscale(self):
        zmin = na.nanmin(self._axes.images[-1]._A)
        zmax = na.nanmax(self._axes.images[-1]._A)
        self.set_zlim(zmin, zmax)

    def switch_y(self, *args, **kwargs):
        pass

    def switch_x(self, *args, **kwargs):
        pass

    def switch_z(self, field):
        self.set_log_field(self.pf.field_info[field].take_log)
        self.axis_names["Z"] = field
        self._redraw_image()

    def selfSetup(self):
        pass

    def _autoset_label(self):
        if self.datalabel is None:
            field_name = self.axis_names["Z"]
            proj = "Proj" in self._type_name and \
                   self.data._weight is None
            data_label = self.pf.field_info[field_name].get_label(proj)
        else: data_label = self.datalabel
        if self.colorbar != None:
            self.colorbar.set_label(str(data_label), **self.label_kws)


class FixedResolutionPlot(VMPlot):

    # This is a great argument in favor of changing the name
    # from VMPlot to something else

    _type_name = "FixedResolution"

    def _get_buff(self, width=None):
        return self.data[self.axis_names["Z"]]

    def set_width(self, width, unit):
        #mylog.debug("Not changing FixedResolution width")
        self._refresh_display_width()

    def _refresh_display_width(self, width=None):
        self.set_xlim(self.data.bounds[0], self.data.bounds[1])
        self.set_ylim(self.data.bounds[2], self.data.bounds[3])
        self._redraw_image()

    def setup_domain_edges(self, *args, **kwargs):
        return

class PCSlicePlot(VMPlot):
    _type_name = "Slice"

    def show_velocity(self, factor = 16, bv_radius = None):
        xax = x_dict[self.data.axis]
        yax = y_dict[self.data.axis]
        xf = "%s-velocity" % (axis_names[xax])
        yf = "%s-velocity" % (axis_names[yax])
        if bv_radius is not None:
            sp = self.data.pf.h.sphere(self.data.center, bv_radius)
            bv = sp.quantities["BulkVelocity"]()
            self.data[xf + "_bv"] = self.data[xf] - bv[xax]
            self.data[yf + "_bv"] = self.data[yf] - bv[yax]
            xf += "_bv"; yf += "_bv"
        from Callbacks import QuiverCallback
        self.add_callback(QuiverCallback(xf, yf, factor))

class NNVMPlot:
    def _get_buff(self, width=None):
        import yt.utilities.delaunay as de
        x0, x1 = self.xlim
        y0, y1 = self.ylim
        if width is None:
            l, b, width, height = mpl_get_bounds(self._axes.bbox)
        else:
            height = width
        self.pix = (width,height)
        numPoints_x = int(width)
        numPoints_y = int(width)
        dx = numPoints_x / (x1-x0)
        dy = numPoints_y / (y1-y0)
        xlim = na.logical_and(self.data["px"]+2.0*self.data['pdx'] >= x0,
                              self.data["px"]-2.0*self.data['pdx'] <= x1)
        ylim = na.logical_and(self.data["py"]+2.0*self.data['pdy'] >= y0,
                              self.data["py"]-2.0*self.data['pdy'] <= y1)
        wI = na.where(na.logical_and(xlim,ylim))
        xi, yi = na.mgrid[0:numPoints_x, 0:numPoints_y]
        x = (self.data["px"][wI]-x0)*dx
        y = (self.data["py"][wI]-y0)*dy
        z = self.data[self.axis_names["Z"]][wI]
        if self.log_field: z=na.log10(z)
        buff = de.Triangulation(x,y).nn_interpolator(z)(xi,yi)
        buff = buff.clip(z.min(), z.max())
        if self.log_field: buff = 10**buff
        return buff.transpose()


class PCSlicePlotNaturalNeighbor(NNVMPlot, PCSlicePlot):
    _type_name = "NNSlice"

class PCProjectionPlot(VMPlot):

    _type_name = "Projection"

    def switch_z(self, field):
        mylog.warning("Choosing not to change the field of a projection instance")

    def _generate_prefix(self, prefix):
        VMPlot._generate_prefix(self, prefix)
        if self.data._weight is not None:
            self.prefix += "_%s" % (self.data._weight)

class PCProjectionPlotNaturalNeighbor(NNVMPlot, PCProjectionPlot):
    _type_name = "NNProj"

class CuttingPlanePlot(PCSlicePlot):

    _type_name = "CuttingPlane"
    def _get_buff(self, width=None):
        px_min, px_max = self.xlim
        py_min, py_max = self.ylim
        if width is None:
            l, b, width, height = mpl_get_bounds(self._axes.bbox)
        else:
            height = width
        self.pix = (width,height)
        indices = na.argsort(self.data['dx'])[::-1]
        buff = _MPL.CPixelize( self.data['x'], self.data['y'], self.data['z'],
                               self.data['px'], self.data['py'],
                               self.data['pdx'], self.data['pdy'], self.data['pdz'],
                               self.data.center, self.data._inv_mat, indices,
                               self.data[self.axis_names['Z']],
                               int(width), int(width),
                               (px_min, px_max, py_min, py_max)).transpose()
        return buff

    def _refresh_display_width(self, width=None):
        
        if width:
            self.width = width
        else:
            width = self.width
        if iterable(width):
            width_x, width_y = width
        else:
            width_x = width
            width_y = width
        l_edge_x = -width_x/2.0
        r_edge_x = +width_x/2.0
        l_edge_y = -width_y/2.0
        r_edge_y = +width_y/2.0
        self.set_xlim(l_edge_x, r_edge_x) # We have no real limits
        self.set_ylim(l_edge_y, r_edge_y) # At some point, perhaps calculate them?
        self._redraw_image()

class ParticlePlot(RavenPlot):

    _type_name = "ParticlePlot"
    def __init__(self, data, axis, width, p_size=1.0, col='k', stride=1.0,
                 figure = None, axes = None):
        kwargs = {}
        if figure is None: kwargs['size'] = (8,8)
        RavenPlot.__init__(self, data, [], figure, axes, **kwargs)
        self._figure.subplots_adjust(hspace=0, wspace=0, bottom=0.0,
                                    top=1.0, left=0.0, right=1.0)
        self.axis = axis
        self.setup_domain_edges(axis)
        self.axis_names['Z'] = col
        self.modify["particles"](width, p_size, col, stride)
        self.set_width(1,'unitary')

    def _redraw_image(self, *args):
        self._axes.clear() # To help out the colorbar
        self._reset_image_parameters()
        self._run_callbacks()

    def _generate_prefix(self, prefix):
        self.prefix = "_".join([prefix, self._type_name, \
            axis_names[self.axis], self.axis_names['Z']])
        self["Field1"] = self.axis_names["Z"]
        self["Field2"] = None
        self["Field3"] = None

    def set_width(self, width, unit):
        self["Unit"] = str(unit)
        self["Width"] = float(width)
        if isinstance(unit, types.StringTypes):
            unit = self.data.pf[str(unit)]
        self.width = width / unit
        self._refresh_display_width()

    def _refresh_display_width(self, width=None):
        if width:
            self.width = width
        else:
            width = self.width
        if iterable(width):
            width_x, width_y = width
        else:
            width_x = width
            width_y = width
        l_edge_x = self.data.center[x_dict[self.axis]] - width_x/2.0
        r_edge_x = self.data.center[x_dict[self.axis]] + width_x/2.0
        l_edge_y = self.data.center[y_dict[self.axis]] - width_y/2.0
        r_edge_y = self.data.center[y_dict[self.axis]] + width_y/2.0
        self.set_xlim(max(l_edge_x,self.xmin), min(r_edge_x,self.xmax))
        self.set_ylim(max(l_edge_y,self.ymin), min(r_edge_y,self.ymax))
        self._redraw_image()

    def _reset_image_parameters(self):
        self._axes.set_xticks(())
        self._axes.set_yticks(())
        self._axes.set_ylabel("")
        self._axes.set_xlabel("")
        l, b, width, height = mpl_get_bounds(self._axes.bbox)
        self._axes.set_xlim(0, width)
        self._axes.set_ylim(0, height)

class ProfilePlot(RavenPlot):
    _x_label = None
    _y_label = None

    def setup_bins(self, field, func=None):
        if field in self.pf.field_info and self.pf.field_info[field].take_log:
            log_field = True
            if func: func('log')
        else:
            log_field = False
            if func: func('linear')
        mylog.debug("Field: %s, log_field: %s", field, log_field)
        return log_field

    def set_x_label(self, label):
        self._axes.set_xlabel(label)
        self._x_label = label

    def set_y_label(self, label):
        self._axes.set_ylabel(label)
        self._y_label = label

    def _autoset_label(self, field, function, axis):
        label = getattr(self, '_%s_label' % axis)
        if label is None: label = self.pf.field_info[field].get_label()
        function(label)

class Profile1DPlot(ProfilePlot):

    _type_name = "Profile1D"
    def __init__(self, data, fields, id, ticker=None, 
                 figure=None, axes=None, plot_options=None):
        self._semi_unique_id = id
        RavenPlot.__init__(self, data, fields, figure, axes)

        self.axis_names["X"] = fields[0]
        self.axis_names["Y"] = fields[1]

        if plot_options is None: plot_options = {}
        self.plot_options = plot_options

        self._log_x = self.data._x_log
        self._log_y = self.setup_bins(self.fields[1])

    def _generate_prefix(self, prefix):
        self.prefix = "_".join([prefix, self._type_name,
                       str(self._semi_unique_id),
                       self.axis_names['X'], self.axis_names['Y']])
        self["Field1"] = self.axis_names["X"]
        self["Field2"] = self.axis_names["Y"]

    def _redraw_image(self):
        self._axes.clear()
        if not self._log_x and not self._log_y:
            func = self._axes.plot
        elif self._log_x and not self._log_y:
            func = self._axes.semilogx
        elif not self._log_x and self._log_y:
            func = self._axes.semilogy
        elif self._log_x and self._log_y:
            func = self._axes.loglog
        indices = na.argsort(self.data[self.fields[0]])
        func(self.data[self.fields[0]][indices],
             self.data[self.fields[1]][indices],
             **self.plot_options)
        self._autoset_label(self.fields[0], self.set_x_label, 'x')
        self._autoset_label(self.fields[1], self.set_y_label, 'y')
        if self.xlim is not None: self._axes.set_xlim(*self.xlim)
        if self.ylim is not None: self._axes.set_ylim(*self.ylim)
        self._run_callbacks()

    def set_log_field(self, val):
        if val:
            self._log_y = True
        else:
            self._log_y = False

    def switch_x(self, field, weight="CellMassMsun", accumulation=False):
        self.fields[0] = field
        self.axis_names["X"] = field
        if field not in self.data.keys():
            self.data.add_fields(field, weight, accumulation)
        self._log_x = self.setup_bins(self.fields[0])
    
    def switch_z(self, field, weight="CellMassMsun", accumulation=False):
        self.fields[1] = field
        self.axis_names["Y"] = field
        if field not in self.data.keys():
            self.data.add_fields(field, weight, accumulation)
        self._log_y = self.setup_bins(self.fields[1])
    switch_y = switch_z # Compatibility...

class PhasePlot(ProfilePlot):

    _type_name = "Profile2D"
    _xlim = None
    _ylim = None
    _z_label = None

    def __init__(self, data, fields, id, ticker=None, cmap=None,
                 figure=None, axes=None):
        self._semi_unique_id = id
        RavenPlot.__init__(self, data, fields, figure, axes)
        self.ticker = ticker
        self.image = None
        self.set_cmap(cmap)
        self._zlim = None
     
        self.axis_names["X"] = fields[0]
        self.axis_names["Y"] = fields[1]
        self.axis_names["Z"] = fields[2]

        self.x_bins = self.data[self.fields[0]]
        self.y_bins = self.data[self.fields[1]]
        self._log_x = self.data._x_log
        self._log_y = self.data._y_log
        self._log_z = self.setup_bins(self.fields[2])
        self.__init_colorbar()

    def _run_callbacks(self):
        # We sublcass to avoid the problem of replacing self.image,
        # which is a collection
        self._axes.patches = []
        self._axes.texts = []
        for cb in self._callbacks:
            cb(self)

    def __init_colorbar(self):
        temparray = na.ones((self.x_bins.size, self.y_bins.size))
        self.norm = matplotlib.colors.Normalize()
        self.image = self._axes.pcolormesh(self.x_bins, self.y_bins,
                                      temparray, shading='flat',
                                      norm=self.norm, cmap=self.cmap,
                                      rasterized=True)
        self.colorbar = self._figure.colorbar(self.image,
                                    extend='neither', shrink=0.95,
                                    format="%0.2e" )

    def set_cmap(self, cmap):
        RavenPlot.set_cmap(self, cmap)
        if self.image != None and self.cmap != None:
            self.image.set_cmap(self.cmap)

    def switch_z(self, field, weight="CellMassMsun", accumulation=False, fractional=False):
        self.fields[2] = field
        self.axis_names["Z"] = field
        if field not in self.data.keys(): self.data.add_fields(field, weight, accumulation, fractional=fractional)
        self._log_z = self.setup_bins(self.fields[2])

    def set_xlim(self, xmin, xmax):
        self._xlim = (xmin,xmax)

    def set_ylim(self, ymin, ymax):
        self._ylim = (ymin,ymax)

    def set_zlim(self, zmin, zmax, dex=None):
        """
        Set the z boundaries of this plot.
        """
        # This next call fixes some things, but is slower...
        #self._redraw_image()
        if (zmin is None) or (zmax is None):    
            if zmin == 'min':
                zmin = na.nanmin(self._axes.images[-1]._A)
                if dex is not None:
                    zmax = min(zmin*10**(dex),na.nanmax(self._axes.images[-1]._A))
            if zmax == 'max':
                zmax = na.nanmax(self._axes.images[-1]._A)
                if dex is not None:
                    zmin = max(zmax/(10**(dex)),na.nanmin(self._axes.images[-1]._A))
        self._zlim = (zmin, zmax)

    def set_log_field(self, val):
        if val:
            self._log_z = True
            self.norm = matplotlib.colors.LogNorm()
            ttype = matplotlib.ticker.LogFormatter
        else:
            self._log_z = False
            self.norm = matplotlib.colors.Normalize()
            ttype = matplotlib.ticker.ScalarFormatter
        if self.colorbar:
            self.colorbar.set_norm(self.norm)
            self.colorbar.formatter = ttype()

    def _redraw_image(self):
        vals = self.data[self.fields[2]].transpose()
        used_bin = self.data["UsedBins"].transpose()
        vmin = na.nanmin(vals[used_bin])
        vmax = na.nanmax(vals[used_bin])
        if self._zlim is not None: vmin, vmax = self._zlim
        if self._log_z:
            # We want smallest non-zero vmin
            self.norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax,
                                                clip=False)
            self.ticker = matplotlib.ticker.LogLocator()
            if self._zlim is None:
                vI = na.where(vals > 0)
                vmin = vals[vI].min()
                vmax = vals[vI].max()
            self.norm.autoscale(na.array((vmin,vmax), dtype='float64'))
        else:
            self.norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax,
                                                  clip=False)
            self.ticker = matplotlib.ticker.MaxNLocator()
        self.colorbar.set_norm(self.norm)
        if self.cmap == None:
            self.cmap = matplotlib.cm.get_cmap()
        self.cmap.set_bad("w")
        self.cmap.set_under("w")
        self.cmap.set_over("w")
        self._axes.clear()
        self._axes.set_xscale("linear")
        self._axes.set_yscale("linear")
        self.image = self._axes.pcolormesh(self.x_bins, self.y_bins, \
                                      vals, shading='flat', \
                                      norm=self.norm, cmap=self.cmap,
                                      rasterized=True)
        self._axes.set_xscale({0:"linear",1:"log"}[int(self._log_x)])
        self._axes.set_yscale({0:"linear",1:"log"}[int(self._log_y)])
        if self._xlim is not None: self._axes.set_xlim(*self._xlim)
        if self._ylim is not None: self._axes.set_ylim(*self._ylim)
        self.vals = vals

        mpl_notify(self.image, self.colorbar)
        self._autoset_label(self.fields[0], self.set_x_label, 'x')
        self._autoset_label(self.fields[1], self.set_y_label, 'y')
        self._autoset_label(self.fields[2], self.set_z_label, 'z')
        self._run_callbacks()

    def _generate_prefix(self, prefix):
        self.prefix = "_".join([str(prefix), self._type_name,
            str(self._semi_unique_id),
            self.axis_names['X'], self.axis_names['Y'],
            self.axis_names['Z'], ])
        self["Field1"] = self.axis_names["X"]
        self["Field2"] = self.axis_names["Y"]
        self["Field3"] = self.axis_names["Z"]

    def set_width(self, width, unit):
        mylog.warning("Choosing not to change the width of a phase plot instance")

    def set_z_label(self, label):
        self.colorbar.set_label(label)
        self._z_label = label

class LineQueryPlot(RavenPlot):
    _type_name = "LineQueryPlot"

    def __init__(self, data, fields, id, ticker=None,
                 figure=None, axes=None, plot_options=None):
        self._semi_unique_id = id
        RavenPlot.__init__(self, data, fields, figure, axes)

        self.axis_names["X"] = fields[0]
        self.axis_names["Y"] = fields[1]

        self._log_x = False
        if fields[1] in self.pf.field_info and \
            self.pf.field_info[fields[1]].take_log:
            self._log_y = True
        else:
            self._log_y = False

        if plot_options is None: plot_options = {}
        self.plot_options = plot_options

    def _generate_prefix(self, prefix):
        self.prefix = "_".join([prefix, self._type_name,
                       str(self._semi_unique_id),
                       self.axis_names['X'], self.axis_names['Y']])
        self["Field1"] = self.axis_names["X"]
        self["Field2"] = self.axis_names["Y"]

    def _redraw_image(self):
        self._axes.clear()
        if not self._log_x and not self._log_y:
            func = self._axes.plot
        elif self._log_x and not self._log_y:
            func = self._axes.semilogx
        elif not self._log_x and self._log_y:
            func = self._axes.semilogy
        elif self._log_x and self._log_y:
            func = self._axes.loglog
        indices = na.argsort(self.data[self.fields[0]])
        func(self.data[self.fields[0]][indices],
             self.data[self.fields[1]][indices],
             **self.plot_options)
        self._autoset_label(self.fields[0], self._axes.set_xlabel)
        self._autoset_label(self.fields[1], self._axes.set_ylabel)
        self._run_callbacks()

    def set_log_field(self, val):
        if val:
            self._log_y = True
        else:
            self._log_y = False

    def switch_x(self, field):
        self.fields[0] = field
        self.axis_names["X"] = field
    
    def switch_z(self, field):
        self.fields[1] = field
        self.axis_names["Y"] = field

    def _autoset_label(self, field_name, func):
        if self.datalabel != None:
            func(str(self.datalabel))
            return
        data_label = r"$\rm{%s}" % field_name.replace("_"," ")
        if field_name in self.pf.field_info:
            data_label += r"\/\/ (%s)" % (self.pf.field_info[field_name].get_units())
        data_label += r"$"
        func(str(data_label))


    switch_y = switch_z # Compatibility...

class ScatterPlot(LineQueryPlot):
    _type_name = "ScatterPlot"

    def _redraw_image(self):
        self._axes.clear()
        self._axes.scatter(
             self.data[self.fields[0]], self.data[self.fields[1]],
             **self.plot_options)
        xscale = {True: "log", False: "linear"}.get(self._log_x)
        yscale = {True: "log", False: "linear"}.get(self._log_y)
        self._axes.set_xscale(xscale)
        self._axes.set_yscale(yscale)
        self._autoset_label(self.fields[0], self._axes.set_xlabel)
        self._autoset_label(self.fields[1], self._axes.set_ylabel)
        self._run_callbacks()

    # We override because we don't want to ditch our collections
    def _run_callbacks(self):
        self._axes.patches = []
        self._axes.texts = []
        for cb in self._callbacks:
            cb(self)

