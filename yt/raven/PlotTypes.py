"""
This is an interface to U{MatPlotLib <http://matplotlib.sf.net>} to plot
irregularly shaped grids, with the presumption that at any point we could have
data that is "hidden" in deeper levels of refinement.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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

from yt.raven import *
from yt.funcs import *

# We only get imported if matplotlib was imported successfully

import _MPL

import matplotlib.image
import matplotlib.ticker
import matplotlib.axes
import matplotlib.figure
import matplotlib._image
import matplotlib.colors
import matplotlib.colorbar
import matplotlib.cm
import matplotlib.collections

def ClusterFilePlot(cls, x, y, xlog=None, ylog=None, fig=None, filename=None,
                    format="png", xbounds = None, ybounds = None):
    """

    """
    if not fig:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        fig = matplotlib.figure.Figure(figsize=(8,8))
        canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot(111)
    if not iterable(cls):
        cls = [cls]
    if xlog == None:
        if lagos.CFfieldInfo.has_key(x):
            xlog = lagos.CFfieldInfo[x][2]
    if ylog == None:
        if lagos.CFfieldInfo.has_key(y):
            ylog = lagos.CFfieldInfo[y][2]
    if xlog and ylog:
        pp=ax.loglog
    elif xlog and not ylog:
        pp=ax.semilogx
    elif ylog and not xlog:
        pp=ax.semilogy
    else:
        pp=ax.plot

    fig.hold(True)
    colors = 'krbgm' * 10
    for cl, cc in zip(cls, colors):
        #pp(cl[x],cl[y], lw=2.5)
        pp(cl[x], cl[y], lw=2.5, color=cc)
    if lagos.CFfieldInfo.has_key(x):
        ax.set_xlabel(lagos.CFfieldInfo[x][1], fontsize=18)
        print lagos.CFfieldInfo[x][1]
    if lagos.CFfieldInfo.has_key(y):
        ax.set_ylabel(lagos.CFfieldInfo[y][1], fontsize=18)
        print lagos.CFfieldInfo[y][1]
    if xbounds:
        ax.set_xlim(xbounds)
    if ybounds:
        ax.set_ylim(ybounds)
    ax.axesFrame.set_linewidth(2)
    for tickLabel in ax.get_xticklabels() + ax.get_yticklabels():
        tickLabel.set_fontsize(14)
    if filename:
        canvas.print_figure(filename, format=format)
    return fig

engineVals = {}

def Initialize(*args, **kwargs):
    engineVals["initialized"] = True
    if not kwargs.has_key("canvas"):
        from matplotlib.backends.backend_agg \
                import FigureCanvasAgg as FigureCanvas
    else:
        FigureCanvas = kwargs["canvas"]
    engineVals["canvas"] = FigureCanvas
    return

def CleanUp(*args, **kwargs):
    pass

class RavenPlot:
    def __init__(self, data, fields, figure = None, axes=None, size=(10,8)):
        self.data = data
        self.fields = fields
        self.size = size
        self.set_autoscale(True)
        self.im = defaultdict(lambda: "")
        self["ParameterFile"] = "%s" % self.data.pf
        self.axis_names = {}
        self._ax_max = self.data.pf["DomainRightEdge"]
        if not figure:
            self._figure = matplotlib.figure.Figure(size)
        else:
            self._figure = figure
        if not figure:
            self._axes = self._figure.add_subplot(1,1,1)
        else:
            self._axes = axes
        self._callbacks = []

    def set_autoscale(self, val):
        self.do_autoscale = val

    def __getitem__(self, item):
        return self.data[item] # Should be returned in CGS

    def save_image(self, prefix, format, submit=None, override=False):
        """
        Save this plot image.  Will generate a filename based on the prefix,
        format, and the approriate data stored in the plot.

        @param prefix: the prefix to prepend to the filename
        @type prefix: string
        @param format: the prefix to append to the filename
        @type format: string
        """
        self._redraw_image()
        if not override:
            self._generate_prefix(prefix)
            my_prefix = self.prefix
        else:
            my_prefix = prefix
        fn = ".".join([my_prefix, format])
        canvas = engineVals["canvas"](self._figure)
        #self._figure.savefig(fn, format)
        canvas.print_figure(fn)
        self["Type"] = self._type_name
        self["GeneratedAt"] = self.data.hierarchy["CurrentTimeIdentifier"]
        return fn

    def _redraw_image(self):
        pass

    def _generate_prefix(self):
        pass

    def set_xlim(self, xmin, xmax):
        self._axes.set_xlim(xmin, xmax)

    def set_ylim(self, ymin, ymax):
        self._axes.set_ylim(ymin, ymax)

    def set_zlim(self, zmin, zmax):
        self._axes.set_zlim(zmin, zmax)

    def set_cmap(self, cmap):
        if isinstance(cmap, types.StringTypes):
            if str(cmap) in raven_colormaps:
                cmap = raven_colormaps[str(cmap)]
            elif hasattr(matplotlib.cm, cmap):
                cmap = getattr(matplotlib.cm, cmap)
        self.cmap = cmap

    def __setitem__(self, item, val):
        self.im[item] = val

    def add_callback(self, func):
        self._callbacks.append(func)
        return len(self._callbacks)-1

    def remove_callback(self, id):
        self._callbacks[id] = lambda a: None

    def _run_callbacks(self):
        for cb in self._callbacks:
            cb(self)

class VMPlot(RavenPlot):
    datalabel = None
    colorbar = None
    def __init__(self, data, field, figure = None, axes = None,
                 use_colorbar = True, size=None):
        fields = ['X', 'Y', field, 'X width', 'Y width']
        if not size:
            size = (10,8)
            if not use_colorbar: size=(8,8)
        RavenPlot.__init__(self, data, fields, figure, axes, size=size)
        self._figure.subplots_adjust(hspace=0, wspace=0, bottom=0.0,
                                    top=1.0, left=0.0, right=1.0)
        self.xmin = 0.0
        self.ymin = 0.0
        self.xmax = 1.0
        self.ymax = 1.0
        self.cmap = None
        if self.data.axis < 3:
            self._x_max = self._ax_max[lagos.x_dict[self.data.axis]]
            self._y_max = self._ax_max[lagos.y_dict[self.data.axis]]
        self.__setup_from_field(field)
        self.__init_temp_image(use_colorbar)

    def __setup_from_field(self, field):
        self.set_log_field(field in lagos.log_fields
                           or lagos.fieldInfo[field].take_log)
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
        else:
            self.colorbar = None
        self.set_width(1,'1')

    def _get_buff(self):
        x0, x1 = self.xlim
        y0, y1 = self.ylim
        l, b, width, height = self._axes.bbox.get_bounds()
        self.pix = (width,height)
        # 'px' == pixel x, or x in the plane of the slice
        # 'x' == actual x
        buff = _MPL.Pixelize(self.data['px'],
                            self.data['py'],
                            self.data['pdx'],
                            self.data['pdy'],
                            self[self.axis_names["Z"]],
                            int(width), int(width),
                        (x0, x1, y0, y1),).transpose()
        return buff

    def _redraw_image(self, *args):
        self._axes.clear() # To help out the colorbar
        buff = self._get_buff()
        mylog.debug("Received buffer of min %s and max %s (%s %s)",
                    buff.min(), buff.max(),
                    self[self.axis_names["Z"]].min(),
                    self[self.axis_names["Z"]].max())
        if self.log_field:
            bI = na.where(buff > 0)
            newmin = buff[bI].min()
            newmax = buff[bI].max()
        else:
            newmin = buff.min()
            newmax = buff.max()
        if self.do_autoscale:
            self.norm.autoscale(na.array((newmin,newmax)))
        self.image = \
            self._axes.imshow(buff, interpolation='nearest', norm = self.norm,
                            aspect=1.0, picker=True, origin='lower')
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
            if self.do_autoscale: self.colorbar.notify(self.image)
        self.autoset_label()

    def set_xlim(self, xmin, xmax):
        self.xlim = (xmin,xmax)

    def set_ylim(self, ymin, ymax):
        self.ylim = (ymin,ymax)

    def _generate_prefix(self, prefix):
        self.prefix = "_".join([prefix, self._type_name, \
            lagos.axis_names[self.data.axis], self.axis_names['Z']])
        self["Field1"] = self.axis_names["Z"]
        self["Field2"] = None
        self["Field3"] = None

    def set_width(self, width, unit):
        self["Unit"] = str(unit)
        self["Width"] = float(width)
        if isinstance(unit, types.StringType):
            unit = self.data.hierarchy[unit]
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
        l_edge_x = self.data.center[lagos.x_dict[self.data.axis]] - width_x/2.0
        r_edge_x = self.data.center[lagos.x_dict[self.data.axis]] + width_x/2.0
        l_edge_y = self.data.center[lagos.y_dict[self.data.axis]] - width_y/2.0
        r_edge_y = self.data.center[lagos.y_dict[self.data.axis]] + width_y/2.0
        self.set_xlim(max(l_edge_x,0.0), min(r_edge_x,self._x_max))
        self.set_ylim(max(l_edge_y,0.0), min(r_edge_y,self._y_max))
        self._redraw_image()

    def autoscale(self):
        zmin = self._axes.images[-1]._A.min()
        zmax = self._axes.images[-1]._A.max()
        self.set_zlim(zmin, zmax)

    def switch_y(self, *args, **kwargs):
        pass

    def switch_x(self, *args, **kwargs):
        pass

    def switch_z(self, field):
        self.set_log_field(field in lagos.log_fields
                           or lagos.fieldInfo[field].take_log)
        self.axis_names["Z"] = field
        self._redraw_image()

    def set_zlim(self, zmin, zmax):
        self.norm.autoscale(na.array([zmin,zmax]))
        self.image.changed()
        if self.colorbar != None:
            self.colorbar.notify(self.image)

    def set_label(self, label):
        if self.colorbar != None: self.colorbar.set_label(str(label))

    def selfSetup(self):
        pass

class SlicePlot(VMPlot):
    _type_name = "Slice"

    def autoset_label(self):
        if self.datalabel != None:
            self.colorbar.set_label(str(self.datalabel))
            return
        field_name = self.axis_names["Z"]
        data_label = r"$\rm{%s}" % field_name
        if lagos.fieldInfo.has_key(field_name):
            data_label += r"\/\/ (%s)" % (lagos.fieldInfo[field_name].get_units())
        data_label += r"$"
        if self.colorbar != None: self.colorbar.set_label(str(data_label))

class ProjectionPlot(VMPlot):
    _type_name = "Projection"

    def autoset_label(self):
        if self.datalabel != None:
            self.colorbar.set_label(str(self.datalabel))
            return
        field_name = self.axis_names["Z"]
        data_label = r"$\rm{%s}" % field_name
        if lagos.fieldInfo.has_key(field_name):
            data_label += r"\/\/ (%s)" % (lagos.fieldInfo[field_name].get_projected_units())
        data_label += r"$"
        if self.colorbar != None: self.colorbar.set_label(str(data_label))

class CuttingPlanePlot(SlicePlot):
    _type_name = "CuttingPlane"
    def _get_buff(self):
        px_min, px_max = self.xlim
        py_min, py_max = self.ylim
        l, b, width, height = self._axes.bbox.get_bounds()
        self.pix = (width,height)
        indices = na.argsort(self.data['dx'])[::-1]
        buff = _MPL.CPixelize( self.data['x'], self.data['y'], self.data['z'],
                               self.data['px'], self.data['py'],
                               self.data['pdx'], self.data['pdy'], self.data['pdz'],
                               self.data.center, self.data._inv_mat, indices,
                               self.data[self.axis_names['Z']],
                               int(width), int(width),
                               (px_min, px_max, py_min, py_max))
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

class PhasePlot(RavenPlot):
    def __init__(self, data, fields, id, ticker=None, cmap=None,
                 figure=None, axes=None):
        self._type_name = "Phase"
        self._semi_unique_id = id
        RavenPlot.__init__(self, data, fields, figure, axes)
        self.ticker = ticker
        self.image = None
        self.set_cmap(cmap)

        self.axis_names["X"] = fields[0]
        self.axis_names["Y"] = fields[1]
        self.axis_names["Z"] = fields[2]

        self.x_bins = self.data[self.fields[0]]
        self.y_bins = self.data[self.fields[1]]
        self._log_x = self.data._x_log
        self._log_y = self.data._y_log
        self._log_z = self.setup_bins(self.fields[2])
        self.__init_colorbar()

    def __init_colorbar(self):
        temparray = na.ones((self.x_bins.size, self.y_bins.size))
        self.norm = matplotlib.colors.Normalize()
        self.image = self._axes.pcolormesh(self.x_bins, self.y_bins,
                                      temparray, shading='flat',
                                      norm=self.norm)
        self.colorbar = self._figure.colorbar(self.image,
                                    extend='neither', shrink=0.95,
                                    format="%0.2e" )

    def setup_bins(self, field, func=None):
        if field in lagos.fieldInfo and lagos.fieldInfo[field].take_log:
            log_field = True
            if func: func('log')
        else:
            log_field = False
            if func: func('linear')
        mylog.debug("Field: %s, log_field: %s", field, log_field)
        return log_field

    def autoset_label(self, field, func):
        dataLabel = r"$\rm{%s}" % (field)
        if field in lagos.fieldInfo:
            dataLabel += r" (%s)" % (lagos.fieldInfo[field].get_units())
        dataLabel += r"$"
        func(str(dataLabel))

    def set_cmap(self, cmap):
        RavenPlot.set_cmap(self, cmap)
        if self.image != None and self.cmap != None:
            self.image.set_cmap(self.cmap)

    def switch_z(self, field, weight="CellMassMsun", accumulation=False):
        self.fields[2] = field
        self.axis_names["Z"] = field
        if field not in self.data.keys(): self.data.add_fields(field, weight, accumulation)
        self._log_z = self.setup_bins(self.fields[2])

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
        vals[~used_bin] = na.nan
        vmin = na.nanmin(vals[used_bin])
        vmax = na.nanmax(vals[used_bin])
        if self._log_z:
            # We want smallest non-zero vmin
            self.norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax,
                                                clip=False)
            self.ticker = matplotlib.ticker.LogLocator()
            vI = na.where(vals > 0)
            newmin = vals[vI].min()
            newmax = vals[vI].max()
            self.norm.autoscale(na.array((newmin,newmax)))
        else:
            self.norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax,
                                                  clip=False)
            self.ticker = matplotlib.ticker.MaxNLocator()
        self.colorbar.set_norm(self.norm)
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
                                      norm=self.norm, cmap=self.cmap)
        self._axes.set_xscale({0:"linear",1:"log"}[int(self._log_x)])
        self._axes.set_yscale({0:"linear",1:"log"}[int(self._log_y)])
        self.vals = vals

        self.colorbar.notify(self.image)
        self.autoset_label(self.fields[0], self._axes.set_xlabel)
        self.autoset_label(self.fields[1], self._axes.set_ylabel)
        self.autoset_label(self.fields[2], self.colorbar.set_label)

    def _generate_prefix(self, prefix):
        self.prefix = "_".join([prefix, self._type_name,
            str(self._semi_unique_id),
            self.axis_names['X'], self.axis_names['Y'],
            self.axis_names['Z'], ])
        self["Field1"] = self.axis_names["X"]
        self["Field2"] = self.axis_names["Y"]
        self["Field3"] = self.axis_names["Z"]

def quiverCallback(field_x, field_y, axis, factor):
    def runCallback(plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        plot._axes.hold(True)
        numPoints_x = plot.image._A.shape[0] / factor
        numPoints_y = plot.image._A.shape[1] / factor
        pixX = _MPL.Pixelize(plot.data['px'],
                             plot.data['py'],
                             plot.data['pdx'],
                             plot.data['pdy'],
                             plot.data[field_x],
                             int(numPoints_x), int(numPoints_y),
                           (x0, x1, y0, y1),).transpose()
        pixY = _MPL.Pixelize(plot.data['px'],
                             plot.data['py'],
                             plot.data['pdx'],
                             plot.data['pdy'],
                             plot.data[field_y],
                             int(numPoints_x), int(numPoints_y),
                           (x0, x1, y0, y1),).transpose()
        X = na.mgrid[0:plot.image._A.shape[0]-1:numPoints_x*1j]# + 0.5*factor
        Y = na.mgrid[0:plot.image._A.shape[1]-1:numPoints_y*1j]# + 0.5*factor
        plot._axes.quiver(X,Y, pixX, -pixY)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)
    return runCallback

def particleCallback(axis, width, p_size=1.0, col='k'):
    field_x = "particle_position_%s" % lagos.axis_names[lagos.x_dict[axis]]
    field_y = "particle_position_%s" % lagos.axis_names[lagos.y_dict[axis]]
    field_z = "particle_position_%s" % lagos.axis_names[axis]
    def runCallback(plot):
        z0 = plot.data.center[axis] - width/2.0
        z1 = plot.data.center[axis] + width/2.0
        grids = plot.data._grids
        particles_x = na.concatenate([g[field_x] for g in grids]).ravel()
        particles_y = na.concatenate([g[field_y] for g in grids]).ravel()
        particles_z = na.concatenate([g[field_z] for g in grids]).ravel()
        if len(particles_x) == 0: return
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        # Now we rescale because our axes limits != data limits
        goodI = na.where( (particles_x < x1) & (particles_x > x0)
                        & (particles_y < y1) & (particles_y > y0)
                        & (particles_z < z1) & (particles_z > z0))
        particles_x = (particles_x[goodI] - x0) * (xx1-xx0)/(x1-x0)
        particles_y = (particles_y[goodI] - y0) * (yy1-yy0)/(y1-y0)
        plot._axes.hold(True)
        plot._axes.scatter(particles_x, particles_y, edgecolors='None',
                          s=p_size, c=col)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)
    return runCallback

def contourCallback(field, ncont=5, factor=4, take_log=False, clim=None):
    try:
        import delaunay as de
    except ImportError:
        mylog.warning("Callback failed; no delaunay module")
        return lambda a: None
    if take_log and clim is not None: clim = (na.log10(clim[0]), na.log10(clim[1]))
    if clim is not None: ncont = na.linspace(clim[0], clim[1], ncont)
    def runCallback(plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        plot._axes.hold(True)
        numPoints_x = plot.image._A.shape[0]
        numPoints_y = plot.image._A.shape[1]
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        xlim = na.logical_and(plot.data["px"] >= x0*0.9,
                              plot.data["px"] <= x1*1.1)
        ylim = na.logical_and(plot.data["py"] >= y0*0.9,
                              plot.data["py"] <= y1*1.1)
        wI = na.where(na.logical_and(xlim,ylim))
        xi, yi = na.mgrid[0:numPoints_x:numPoints_x/(factor*1j),\
                          0:numPoints_y:numPoints_y/(factor*1j)]
        x = (plot.data["px"][wI]-x0)*dx
        y = (plot.data["py"][wI]-y0)*dy
        z = plot.data[field][wI]
        if take_log: z=na.log10(z)
        zi = de.Triangulation(x,y).nn_interpolator(z)(xi,yi)
        print ncont, zi.min(), zi.max()
        plot._axes.contour(xi,yi,zi,ncont,colors='k')
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)
    return runCallback

def gridBoundaryCallback(alpha=1.0, min_pix = 1):
    def runCallback(plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        GLE = plot.data.gridLeftEdge
        GRE = plot.data.gridRightEdge
        px_index = lagos.x_dict[plot.data.axis]
        py_index = lagos.y_dict[plot.data.axis]
        left_edge_px = (GLE[:,px_index]-x0)*dx
        left_edge_py = (GLE[:,py_index]-y0)*dy
        right_edge_px = (GRE[:,px_index]-x0)*dx
        right_edge_py = (GRE[:,py_index]-y0)*dy
        verts = na.array(
                [(left_edge_px, left_edge_px, right_edge_px, right_edge_px),
                 (left_edge_py, right_edge_py, right_edge_py, left_edge_py)])
        visible =  ( right_edge_px - left_edge_px > min_pix ) & \
                   ( right_edge_px - left_edge_px > min_pix )
        verts=verts.transpose()[visible,:,:]
        edgecolors = (0.0,0.0,0.0,alpha)
        grid_collection = matplotlib.collections.PolyCollection(
                verts, facecolors=(0.0,0.0,0.0,0.0),
                       edgecolors=edgecolors)
        plot._axes.hold(True)
        plot._axes.add_collection(grid_collection)
        plot._axes.hold(False)
    return runCallback

def get_smallest_appropriate_unit(v, pf):
    max_nu = 1e30
    good_u = None
    for unit in ['mpc','kpc','pc','au','rsun','cm']:
        vv = v*pf[unit]
        if vv < max_nu and vv > 1.0:
            good_u = unit
            max_nu = v*pf[unit]
    return good_u

def unitBoundaryCallback(unit = "au", factor=4, text_annotate=True, text_which=-2):
    def runCallback(plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        l, b, width, height = plot._axes.bbox.get_bounds()
        xi = lagos.x_dict[plot.data.axis]
        yi = lagos.y_dict[plot.data.axis]
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        center = plot.data.center
        min_dx = plot.data['pdx'].min()
        max_dx = plot.data['pdx'].max()
        w_min_x = 250.0 * min_dx
        w_max_x = 1.0 / factor
        min_exp_x = na.ceil(na.log10(w_min_x*plot.data.pf[unit])/na.log10(factor))
        max_exp_x = na.floor(na.log10(w_max_x*plot.data.pf[unit])/na.log10(factor))
        n_x = max_exp_x - min_exp_x + 1
        widths = na.logspace(min_exp_x, max_exp_x, num = n_x, base=factor)
        widths /= plot.data.pf[unit]
        left_edge_px = (center[xi] - widths/2.0 - x0)*dx
        left_edge_py = (center[yi] - widths/2.0 - y0)*dy
        right_edge_px = (center[xi] + widths/2.0 - x0)*dx
        right_edge_py = (center[yi] + widths/2.0 - y0)*dy
        verts = na.array(
                [(left_edge_px, left_edge_px, right_edge_px, right_edge_px),
                 (left_edge_py, right_edge_py, right_edge_py, left_edge_py)])
        visible =  ( right_edge_px - left_edge_px > 25 ) & \
                   ( right_edge_px - left_edge_px > 25 ) & \
                   ( (right_edge_px < width) & (left_edge_px > 0) ) & \
                   ( (right_edge_py < height) & (left_edge_py > 0) )
        verts=verts.transpose()[visible,:,:]
        grid_collection = matplotlib.collections.PolyCollection(
                verts, facecolors=(0.0,0.0,0.0,0.0),
                       edgecolors = (0.0,0.0,0.0,1.0),
                       linewidths=2.5)
        plot._axes.hold(True)
        plot._axes.add_collection(grid_collection)
        if text_annotate:
            ti = max(text_which, -1*len(widths[visible]))
            w = widths[visible][ti]
            good_u = get_smallest_appropriate_unit(w, plot.data.pf)
            w *= plot.data.pf[good_u]
            plot._axes.annotate("%0.3e %s" % (w,good_u), verts[ti,1,:]+5)
        plot._axes.hold(False)
    return runCallback
