"""
This is an interface to U{MatPlotLib <http://matplotlib.sf.net>} to plot
irregularly shaped grids, with the presumption that at any point we could have
data that is "hidden" in deeper levels of refinement.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from yt.raven import *
from yt.funcs import *

# We only get imported if matplotlib was imported successfully

import _MPL

import matplotlib.image
import matplotlib.axes
import matplotlib.figure
import matplotlib._image
import matplotlib.colors
import matplotlib.colorbar
import matplotlib.cm
from matplotlib.backends.backend_agg import FigureCanvasAgg

class LinkedAMRSubPlots:
    def __init__(self, linkWidth, linkZ, plots = []):
        pass
    def SaveFigure(self, filename, format):
        pass
    def setWidth(self, width, unit):
        pass
    def setZ(self, zmin, zmax):
        pass
    def __getitem__(self, item):
        pass

class ClassicThreePane:
    pass

def ClusterFilePlot(cls, x, y, xlog=None, ylog=None, fig=None, filename=None,
                    format="png", xbounds = None, ybounds = None):
    """
    
    """
    if not fig:
        fig = matplotlib.figure.Figure(figsize=(8,8))
        canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(hspace=0,wspace=0,bottom=0.0, top=1.0, left=0.0, right=1.0)
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
    for cl in cls:
        pp(cl[x],cl[y])
    if lagos.CFfieldInfo.has_key(x):
        ax.set_xlabel(lagos.CFfieldInfo[x][1])
    if lagos.CFfieldInfo.has_key(y):
        ax.set_ylabel(lagos.CFfieldInfo[y][1])
    if xbounds:
        ax.set_xlim(xbounds)
    if ybounds:
        ax.set_ylim(ybounds)
    if filename:
        canvas.print_figure(filename, format=format)
    else:
        return fig

from collections import defaultdict

engineVals = {}
skipAxes = ["X WIDTH", "Y WIDTH", "WEIGHT (OPTIONAL)", "DX", "DY"]

axisFieldDict = {'X':'Field1', 'Y':'Field2', 'Z':'Field3'}

def Initialize(*args, **kwargs):
    engineVals["initialized"] = True
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    engineVals["canvas"] = FigureCanvas
    return

def CleanUp(*args, **kwargs):
    pass

class RavenPlot:
    def __init__(self, data, fields):
        self.data = data
        self.fields = fields
        # With matplotlib, we don't need to copy data back and forth.
        # Unfortunately, it is also harder to plot things, as they may not have
        # a consistent interface.  Also, the plots are designed to be much more
        # *static* than in hippodraw, so switching axes is going to be a bit of
        # a pain.
        self.im = defaultdict(lambda: "")
        self["ParameterFile"] = \
            self.data.hierarchy.parameterFile.parameterFilename
        self.axisNames = {}
        self.figure = matplotlib.figure.Figure((10,8))
        #self.axes = self.figure.add_axes(aspect='equal')
        self.axes = self.figure.add_subplot(1,1,1)#,aspect=1.0)
        #self.axes.set_aspect(1.0)

    def __getitem__(self, item):
        return self.data[item] * \
                    self.data.hierarchy.parameterFile.conversionFactors[item]

    def saveImage(self, prefix, format, submit=None):
        """
        Save this plot image.  Will generate a filename based on the prefix,
        format, and the approriate data stored in the plot.

        @param prefix: the prefix to prepend to the filename
        @type prefix: string
        @param format: the prefix to append to the filename
        @type format: string
        """
        self.generatePrefix(prefix)
        fn = ".".join([self.prefix, format])
        canvas = engineVals["canvas"](self.figure)
        #self.figure.savefig(fn, format)
        canvas.print_figure(fn,format=format)
        self["Type"] = self.typeName
        self["GeneratedAt"] = self.data.hierarchy["CurrentTimeIdentifier"]
        return fn

    def set_xlim(self, xmin, xmax):
        self.axes.set_xlim(xmin, xmax)
        
    def set_ylim(self, ymin, ymax):
        self.axes.set_ylim(ymin, ymax)

    def set_zlim(self, zmin, zmax):
        self.axes.set_zlim(zmin, zmax)

    def __setitem__(self, item, val):
        self.im[item] = val

class VMPlot(RavenPlot):

    def __init__(self, data, field):
        fields = ['X', 'Y', field, 'X width', 'Y width']
        RavenPlot.__init__(self, data, fields)
        self.figure.subplots_adjust(hspace=0,wspace=0,bottom=0.0, top=1.0, left=0.0, right=1.0)
        self.xmin = 0.0
        self.ymin = 0.0
        self.xmax = 1.0
        self.ymax = 1.0
        self.norm = matplotlib.colors.LogNorm()
        temparray = na.random.rand(640,640)
        self.image = \
            self.axes.imshow(temparray, interpolation='nearest', norm = self.norm,
                            aspect=1.0)
        self.axisNames["Z"] = field
        self.axes.set_xticks(())
        self.axes.set_yticks(())
        self.axes.set_ylabel("")
        self.axes.set_xlabel("")
        self.colorbar = self.figure.colorbar(self.axes.images[-1], \
                                             extend='neither', \
                                             shrink=0.95)
        self.set_width(1,'1')
        self.redraw_image()

    def redraw_image(self):
        #x0, y0, v_width, v_height = self.axes.viewLim.get_bounds()
        x0, x1 = self.xlim
        y0, y1 = self.ylim
        l, b, width, height = self.axes.bbox.get_bounds()
        buff = _MPL.Pixelize(self.data['x'],
                             self.data['y'],
                             self.data['dx'],
                             self.data['dy'],
                             self[self.axisNames["Z"]],
                             int(height), int(width),
                           (x0, x1, y0, y1),)
        self.image.set_data(buff)
        self.norm.autoscale((buff.min(),buff.max()))
        self.colorbar.notify(self.image)

    def set_xlim(self, xmin, xmax):
        self.xlim = (xmin,xmax)
        
    def set_ylim(self, ymin, ymax):
        self.ylim = (ymin,ymax)

    def generatePrefix(self, prefix):
        self.prefix = "_".join([prefix, self.typeName, \
            lagos.axis_names[self.data.axis], self.axisNames['Z']])
        self["Field1"] = self.axisNames["Z"]
        self["Field2"] = ""
        self["Field3"] = ""

    def set_width(self, width, unit):
        self["Unit"] = str(unit)
        self["Width"] = float(width)
        if isinstance(unit, types.StringType):
            unit = self.data.hierarchy[unit]
        self.width = width / unit
        self.refreshDisplayWidth()

    def checkColormap(self, field):
        cmap = lagos.colormap_dict[field]
        #self.image.set_cmap(cmap)

    def refreshDisplayWidth(self, width=None):
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
        self.redraw_image()

    def autoscale(self):
        zmin = self.axes.images[-1]._A.min()
        zmax = self.axes.images[-1]._A.max()
        self.set_zlim(zmin, zmax)

    def switch_y(self, *args, **kwargs):
        pass

    def switch_x(self, *args, **kwargs):
        pass

    def switch_z(self, field):
        self.axisNames["Z"] = field
        self.redraw_image()

    def set_zlim(self, zmin, zmax):
        self.norm.autoscale(na.array([zmin,zmax]))
        self.image.changed()
        self.colorbar.notify(self.image)

    def set_label(self, label):
        self.colorbar.set_label(label)

    def set_cmap(self, cmap):
        if isinstance(cmap, types.StringType):
            if hasattr(matplotlib.cm, cmap):
                cmap = getattr(matplotlib.cm, cmap)
        self.image.set_cmap(cmap)

class SlicePlot(VMPlot):
    def __init__(self, data, fields):
        self.typeName = "Slice"
        VMPlot.__init__(self, data, fields)

    def autoset_label(self):
        dataLabel = self.axisNames["Z"]
        if lagos.fieldInfo.has_key(self.axisNames["Z"]):
            dataLabel += " (%s)" % (lagos.fieldInfo[self.axisNames["Z"]][0])
        self.colorbar.set_label(dataLabel)

class ProjectionPlot(VMPlot):
    def __init__(self, data, fields):
        self.typeName = "Projection"
        VMPlot.__init__(self, data, fields)

    def autoset_label(self):
        dataLabel = self.axisNames["Z"]
        if lagos.fieldInfo.has_key(self.axisNames["Z"]):
            dataLabel += " (%s)" % (lagos.fieldInfo[self.axisNames["Z"]][1])
        self.colorbar.set_label(dataLabel)

class TwoPhasePlot(RavenPlot):
    def __init__(self, data, fields, width=None, unit=None):
        self.typeName = "TwoPhase"
        RavenPlot.__init__(self, data, fields)
        self.axisNames["X"] = fields[0]
        self.axisNames["Y"] = fields[1]
        x_v = self.data[self.axisNames["X"]]
        y_v = self.data[self.axisNames["Y"]]
        self.axes.set_xscale("log")
        self.axes.set_yscale("log")
        x_bins = na.logspace(na.log10(x_v.min()),na.log10(x_v.max()),num=101)
        y_bins = na.logspace(na.log10(y_v.min()),na.log10(y_v.max()),num=101)
        vals, x, y = na.histogram2d( \
            na.log10(x_v), na.log10(y_v),
            bins = 100, normed=True )
        i = na.where(vals>0)
        vmin = vals[i].min()
        self.norm=matplotlib.colors.LogNorm(vmin=vmin, clip=False)
        self.cmap = matplotlib.cm.get_cmap()
        self.cmap.set_under("k")
        self.image = self.axes.pcolor(x_bins,y_bins, \
                                      vals,shading='flat', \
                                      norm=self.norm)
        self.colorbar = self.figure.colorbar(self.image, \
                                             extend='neither', \
                                             shrink=0.95, cmap=self.cmap)
        self.colorbar.notify(self.image)

    def generatePrefix(self, prefix):
        self.prefix = "_".join([prefix, self.typeName, \
            self.axisNames['X'], self.axisNames['Y']])
        self["Field1"] = self.axisNames["X"]
        self["Field2"] = self.axisNames["Y"]
        self["Field3"] = ""

class ThreePhasePlot(RavenPlot):
    def __init__(self, data, fields, width=None, unit=None):
        nbins = 100
        self.typeName = "ThreePhase"
        RavenPlot.__init__(self, data, fields)
        self.axisNames["X"] = fields[0]
        self.axisNames["Y"] = fields[1]
        self.axisNames["Z"] = fields[2]
        x_v = self.data[self.axisNames["X"]]
        y_v = self.data[self.axisNames["Y"]]
        z_v = self.data[self.axisNames["Z"]]
        weight = self.data["CellMass"]
        x_bins = na.logspace(na.log10(x_v.min()*0.99),na.log10(x_v.max()*1.01),num=nbins)
        y_bins = na.logspace(na.log10(y_v.min()*0.99),na.log10(y_v.max()*1.01),num=nbins)
        x_bins_ids = na.digitize(x_v, x_bins)
        y_bins_ids = na.digitize(y_v, y_bins)
        vals = na.ones((nbins,nbins), dtype=nT.Float64)
        weight_vals = na.zeros((nbins,nbins), dtype=nT.Float64)
        for k in range(len(x_v)):
            j,i = x_bins_ids[k], y_bins_ids[k]
            weight_vals[i,j] += weight[k]
            vals[i,j] += z_v[k]*weight[k]
        vals = vals / weight_vals
        vi = na.where(na.isinf(abs(vals)))
        vals[vi] = na.nan
        self.axes.set_xscale("log")
        self.axes.set_yscale("log")
        vmin = na.nanmin(vals)
        vmax = na.nanmax(vals)
        vals[vi] = 0.0
        self.norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax, clip=False)
        self.cmap = matplotlib.cm.get_cmap()
        self.cmap.set_bad("k")
        self.cmap.set_under("k")
        self.cmap.set_over("k")
        self.image = self.axes.pcolor(x_bins,y_bins, \
                                      vals,shading='flat', \
                                      norm=self.norm)
        self.colorbar = self.figure.colorbar(self.image, \
                                             extend='neither', \
                                             shrink=0.95, cmap=self.cmap)
        self.colorbar.notify(self.image)

    def generatePrefix(self, prefix):
        self.prefix = "_".join([prefix, self.typeName, \
            self.axisNames['X'], self.axisNames['Y'], \
            self.axisNames['Z']])
        self["Field1"] = self.axisNames["X"]
        self["Field2"] = self.axisNames["Y"]
        self["Field3"] = self.axisNames["Z"]

