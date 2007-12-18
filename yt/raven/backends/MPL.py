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

from collections import defaultdict

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
skipAxes = ["X WIDTH", "Y WIDTH", "WEIGHT (OPTIONAL)", "DX", "DY"]

axisFieldDict = {'X':'Field1', 'Y':'Field2', 'Z':'Field3'}

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
        if isinstance(cmap, types.StringType):
            if hasattr(matplotlib.cm, cmap):
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
    def __init__(self, data, field, figure = None, axes = None,
                 use_colorbar = True, size=(800,800)):
        fields = ['X', 'Y', field, 'X width', 'Y width']
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
        self.__setup_from_field(field)
        self.__init_temp_image(use_colorbar)

    def __setup_from_field(self, field):
        if field in lagos.log_fields or lagos.fieldInfo[field].take_log:
            self.log_field = True
            self.norm = matplotlib.colors.LogNorm()
        else:
            self.log_field = False
            self.norm = matplotlib.colors.Normalize()
        self.axis_names["Z"] = field

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

    def _redraw_image(self, *args):
        self._axes.clear() # To help out the colorbar
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
        self._axes.set_xticks(())
        self._axes.set_yticks(())
        self._axes.set_ylabel("")
        self._axes.set_xlabel("")
        if self.cmap:
            self.image.set_cmap(self.cmap)
        if self.colorbar != None:
            self.colorbar.notify(self.image)
        self.autoset_label()
        self._run_callbacks()

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
        l_edge_x = self.data.center[lagos.x_dict[self.data.axis]] - width/2.0
        r_edge_x = self.data.center[lagos.x_dict[self.data.axis]] + width/2.0
        l_edge_y = self.data.center[lagos.y_dict[self.data.axis]] - width/2.0
        r_edge_y = self.data.center[lagos.y_dict[self.data.axis]] + width/2.0
        self.set_xlim(max(l_edge_x,0.0), min(r_edge_x,1.0))
        self.set_ylim(max(l_edge_y,0.0), min(r_edge_y,1.0))
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
        if field in lagos.log_fields or lagos.fieldInfo[field].take_log:
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
        self.axis_names["Z"] = field
        self._redraw_image()

    def set_zlim(self, zmin, zmax):
        self.norm.autoscale(na.array([zmin,zmax]))
        self.image.changed()
        if self.colorbar != None:
            self.colorbar.notify(self.image)

    def set_label(self, label):
        if self.colorbar != None: self.colorbar.set_label(label)

    def selfSetup(self):
        pass

class SlicePlot(VMPlot):
    _type_name = "Slice"

    def autoset_label(self):
        if self.datalabel != None:
            self.colorbar.set_label(self.datalabel)
            return
        field_name = self.axis_names["Z"]
        if lagos.fieldInfo.has_key(field_name):
            field_name += " (%s)" % (lagos.fieldInfo[field_name].get_units())
        if self.colorbar != None: self.colorbar.set_label(field_name)

class ProjectionPlot(VMPlot):
    _type_name = "Projection"

    def autoset_label(self):
        dataLabel = self.axis_names["Z"]
        field_name = self.axis_names["Z"]
        if lagos.fieldInfo.has_key(field_name):
            field_name += " (%s)" % (lagos.fieldInfo[field_name].get_projected_units())
        if self.colorbar != None: self.colorbar.set_label(field_name)

    def __getitem__(self, item):
        return self.data[item] * \
                  self.data.hierarchy.parameter_file.units["cm"]
                  #self.data.hierarchy.parameter_file.conversion_factors[item] * \

class PhasePlot(RavenPlot):
    def __init__(self, data, fields, width=None, unit=None, bins=100,
                 weight=None, ticker=None, cmap=None, figure=None, axes=None):
        self._type_name = "Phase"
        RavenPlot.__init__(self, data, fields, figure, axes)
        self.ticker = ticker
        self.image = None
        self.bins = bins
        self.set_cmap(cmap)
        self.weight = weight

        self.axis_names["X"] = fields[0]
        self.axis_names["Y"] = fields[1]
        self.axis_names["Z"] = fields[2]

        log_field, self.x_v, self.x_bins = self.setup_bins(self.fields[0],
                                                       self._axes.set_xscale)
        log_field, self.y_v, self.y_bins = self.setup_bins(self.fields[1],
                                                       self._axes.set_yscale)
        self.log_z, self.z_v, self.z_bins = self.setup_bins(self.fields[2])

        self.colorbar = None

    def setup_bins(self, field, func = None):
        log_field = False
        v = self.data[field]
        if field in lagos.log_fields or lagos.fieldInfo[field].take_log:
            log_field = True
            bins = na.logspace(na.log10(v.min()*0.99),
                               na.log10(v.max()*1.01),
                               num=self.bins)
            if func: func('log')
        else:
            bins = na.linspace(v.min()*0.99,v.max()*1.01,num=self.bins)
            if func: func('linear')
        return log_field, v, bins

    def autoset_label(self, field, func):
        dataLabel = field
        if lagos.fieldInfo.has_key(field):
            dataLabel += " (%s)" % (lagos.fieldInfo[field].get_units())
        func(dataLabel)

    def set_cmap(self, cmap):
        RavenPlot.set_cmap(self, cmap)
        if self.image != None and self.cmap != None:
            self.image.set_cmap(self.cmap)

    def switch_x(self, field):
        self.fields[0] = field
        self.axis_names["X"] = field
        log_field, self.x_v, self.x_bins = self.setup_bins(self.fields[0],
                                                       self._axes.set_xscale)

    def switch_y(self, field):
        self.fields[1] = field
        self.axis_names["Y"] = field
        log_field, self.y_v, self.y_bins = self.setup_bins(self.fields[1],
                                                       self._axes.set_yscale)

    def switch_z(self, field):
        self.fields[2] = field
        self.axis_names["Z"] = field
        self.log_z, self.z_v, self.z_bins = self.setup_bins(self.fields[2])

    def switch_weight(self, weight):
        if weight == "": weight=None
        self.weight = weight

    def _redraw_image(self):
        l, b, width, height = self._axes.bbox.get_bounds()
        self.pix = (width,height)
        x_bins_ids = na.digitize(self.x_v, self.x_bins)
        y_bins_ids = na.digitize(self.y_v, self.y_bins)

        vals = na.zeros((self.bins,self.bins), dtype='float64')
        weight_vals = na.zeros((self.bins,self.bins), dtype='float64')
        used_bin = na.zeros((self.bins,self.bins), dtype='bool')

        x_ind, y_ind = (x_bins_ids-1,y_bins_ids-1) # To match up with pcolor
                # pcolor expects value to be between i and i+1, digitize gives
                # bin between i-1 and i
        used_bin[y_ind,x_ind] = True
        nx = len(self.x_v)
        if self.weight != None:
            weight = self.data[self.weight]
        else:
            weight = na.ones(nx)

        z_v = self.z_v
        code =r"""
               int i,j;
               for(int n = 0; n < nx ; n++) {
                 //printf("%d\n",n);
                 j = x_bins_ids(n)-1;
                 i = y_bins_ids(n)-1;
                 weight_vals(i,j) += weight(n);
                 vals(i,j) += z_v(n) * weight(n);
               }
               """
        try:
            weave.inline(code, ['nx','x_bins_ids','y_bins_ids',
                                'weight_vals','weight','vals','z_v'],
                         compiler='gcc', type_converters=converters.blitz,
                         auto_downcast = 0)
        except:
            mylog.debug("SciPy weaving did not work; falling back on loop")
            for k in range(nx):
                j,i = x_bins_ids[k]-1, y_bins_ids[k]-1
                weight_vals[i,j] += weight[k]
                vals[i,j] += self.z_v[k]*weight[k]

        vi = na.where(used_bin == False)
        vit = na.where(used_bin == True)
        if self.weight != None: vals = vals / weight_vals

        vmin = na.nanmin(vals[vit])
        vmax = na.nanmax(vals[vit])
        vals[vi] = 0.0
        if self.log_z:
            # We want smallest non-zero vmin
            vmin=vals[vals>0.0].min()
            self.norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax,
                                                clip=False)
            location_of_ticks = na.logspace(vmin*1.1, vmax*0.9, num=6)
            self.ticker = matplotlib.ticker.LogLocator()
        else:
            self.ticker = matplotlib.ticker.MaxNLocator()
            self.norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax,
                                                  clip=False)
        if self.cmap == None:
            self.cmap = matplotlib.cm.get_cmap()
        self.cmap.set_bad("w")
        self.cmap.set_under("w")
        self.cmap.set_over("w")
        self._axes.clear()
        xs = self._axes.get_xscale()
        ys = self._axes.get_yscale()
        self._axes.set_xscale("linear")
        self._axes.set_yscale("linear")
        self.image = self._axes.pcolormesh(self.x_bins, self.y_bins, \
                                      vals,shading='flat', \
                                      norm=self.norm, cmap=self.cmap)
        self.vals = vals
        self._axes.set_xscale(xs)
        self._axes.set_yscale(ys)
        #self.ticker = matplotlib.ticker.LogLocator(subs=[0.25, 0.5, 0.75, 1])

        if self.colorbar == None:
            self.colorbar = self._figure.colorbar(self.image, \
                                                 extend='neither', \
                                                 shrink=0.95, cmap=self.cmap, \
                                   ticks = self.ticker, format="%0.2e" )

        self.colorbar.notify(self.image)

        self.autoset_label(self.fields[0], self._axes.set_xlabel)
        self.autoset_label(self.fields[1], self._axes.set_ylabel)
        self.autoset_label(self.fields[2], self.colorbar.set_label)

    def _generate_prefix(self, prefix):
        self.prefix = "_".join([prefix, self._type_name, \
            self.axis_names['X'], self.axis_names['Y'], \
            self.axis_names['Z']])
        self["Field1"] = self.axis_names["X"]
        self["Field2"] = self.axis_names["Y"]
        self["Field3"] = self.axis_names["Z"]

def quiverCallback(field_x, field_y, axis, factor):
    def runCallback(plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot.axes.get_xlim()
        yy0, yy1 = plot.axes.get_ylim()
        plot.axes.hold(True)
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
        plot.axes.quiver(X,Y, pixX, -pixY)
        plot.axes.set_xlim(xx0,xx1)
        plot.axes.set_ylim(yy0,yy1)
        plot.axes.hold(False)
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
        xx0, xx1 = plot.axes.get_xlim()
        yy0, yy1 = plot.axes.get_ylim()
        # Now we rescale because our axes limits != data limits
        goodI = na.where( (particles_x < x1) & (particles_x > x0)
                        & (particles_y < y1) & (particles_y > y0)
                        & (particles_z < z1) & (particles_z > z0))
        particles_x = (particles_x[goodI] - x0) * (xx1-xx0)/(x1-x0)
        particles_y = (particles_y[goodI] - y0) * (yy1-yy0)/(y1-y0)
        plot.axes.hold(True)
        plot.axes.scatter(particles_x, particles_y, edgecolors='None',
                          s=p_size, c=col)
        plot.axes.set_xlim(xx0,xx1)
        plot.axes.set_ylim(yy0,yy1)
        plot.axes.hold(False)
    return runCallback

def contourCallback(field, axis, ncont=5, factor=4):
    try:
        import scipy.sandbox.delaunay as de
    except ImportError:
        mylog.warning("Callback failed; no delaunay module")
        return lambda a: None
    def runCallback(plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot.axes.get_xlim()
        yy0, yy1 = plot.axes.get_ylim()
        plot.axes.hold(True)
        numPoints_x = plot.image._A.shape[0]
        numPoints_y = plot.image._A.shape[1]
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        xlim = na.logical_and(plot.data["x"] >= x0*0.9,
                              plot.data["x"] <= x1*1.1)
        ylim = na.logical_and(plot.data["y"] >= y0*0.9,
                              plot.data["y"] <= y1*1.1)
        wI = na.where(na.logical_and(xlim,ylim))
        xi, yi = na.mgrid[0:numPoints_x:numPoints_x/(factor*1j),\
                          0:numPoints_y:numPoints_y/(factor*1j)]
        x = (plot.data["x"][wI]-x0)*dx
        y = (plot.data["y"][wI]-y0)*dy
        z = plot.data[field][wI]
        zi = de.Triangulation(x,y).nn_interpolator(z)(xi,yi)
        plot.axes.contour(xi,yi,zi,ncont, colors='k')
        plot.axes.set_xlim(xx0,xx1)
        plot.axes.set_ylim(yy0,yy1)
        plot.axes.hold(False)
    return runCallback
