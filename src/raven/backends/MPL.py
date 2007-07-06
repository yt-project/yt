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
import matplotlib.ticker
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
        temparray = na.ones((800,800))
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
                             int(800), int(800),
                           (x0, x1, y0, y1),)
        self.norm.autoscale(na.array((buff.min(),buff.max())))
        self.image.set_data(buff)
        self.colorbar.notify(self.image)
        self.autoset_label()

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

    def __getitem__(self, item):
        return self.data[item] * \
                    self.data.hierarchy.parameterFile.conversionFactors[item]


class PhasePlot(RavenPlot):
    def __init__(self, data, fields, bins = 100, width=None, unit=None):
        RavenPlot.__init__(self, data, fields)
        self.bins = bins
        self.axisNames["X"] = fields[0]
        self.axisNames["Y"] = fields[1]
        logIt, self.x_v, self.x_bins = self.setup_bins(fields[0], self.axes.set_xscale)
        logIt, self.y_v, self.y_bins = self.setup_bins(fields[1], self.axes.set_yscale)

    def setup_bins(self, field, func):
        logIt = False
        v = na.log10(self.data[field])
        if field in lagos.log_fields or lagos.fieldInfo[field][2]:
            logIt = True
            bins = na.logspace(na.log10(v.min()),na.log10(v.max()),num=self.bins+1)
            func('log')
        else:
            bins = na.linspace(v.min(),v.max(),num=self.bins+1)
            func('linear')
        return logIt, v, bins

    def autoset_label(self, field, func):
        dataLabel = field
        if lagos.fieldInfo.has_key(field):
            dataLabel += " (%s)" % (lagos.fieldInfo[field][0])
        func(dataLabel)

class TwoPhasePlot(PhasePlot):
    def __init__(self, data, fields, bins = 100, width=None, unit=None):
        self.typeName = "TwoPhase"
        PhasePlot.__init__(self, data, fields, bins, width, unit)

        vals, x, y = na.histogram2d( \
            self.x_v, self.y_v, \
            bins = (self.x_bins, self.y_bins), \
            normed=False )
        i = na.where(vals>0)
        vmin = vals[i].min()
        self.norm=matplotlib.colors.LogNorm(vmin=vmin, clip=False)
        self.cmap = matplotlib.cm.get_cmap()
        self.cmap.set_under("k")
        self.image = self.axes.pcolor(self.x_bins,self.y_bins, \
                                      vals.transpose(),shading='flat', \
                                      norm=self.norm)

        self.autoset_label(fields[0], self.axes.set_xlabel)
        self.autoset_label(fields[1], self.axes.set_ylabel)

        self.colorbar = self.figure.colorbar(self.image, \
                                             extend='neither', \
                                             shrink=0.95, cmap=self.cmap)
        self.colorbar.notify(self.image)
        self.colorbar.set_label("Cells per Bin")

    def generatePrefix(self, prefix):
        self.prefix = "_".join([prefix, self.typeName, \
            self.axisNames['X'], self.axisNames['Y']])
        self["Field1"] = self.axisNames["X"]
        self["Field2"] = self.axisNames["Y"]
        self["Field3"] = None

class ThreePhasePlot(PhasePlot):
    def __init__(self, data, fields, width=None, unit=None, bins=100, weight="CellMass"):
        self.typeName = "ThreePhase"
        PhasePlot.__init__(self, data, fields)

        self.axisNames["Z"] = fields[2]
        logIt, self.z_v, self.z_bins = self.setup_bins(fields[2], lambda i: None)

        weight = self.data[weight]
        x_bins_ids = na.digitize(self.x_v, self.x_bins)
        y_bins_ids = na.digitize(self.y_v, self.y_bins)

        vals = na.zeros((self.bins+1,self.bins+1), dtype=nT.Float64) 
        weight_vals = na.zeros((self.bins+1,self.bins+1), dtype=nT.Float64)
        used_bin = na.zeros((self.bins+1,self.bins+1), dtype=nT.Bool)

        for k in range(len(self.x_v)):
            j,i = x_bins_ids[k]-1, y_bins_ids[k]-1
            used_bin[i,j] = True
            weight_vals[i,j] += weight[k]
            vals[i,j] += self.z_v[k]*weight[k]

        vi = na.where(used_bin == False)
        vit = na.where(used_bin == True)
        vals = vals / weight_vals
        self.axes.set_xscale("log")
        self.axes.set_yscale("log")
        vmin = na.nanmin(vals[vit])
        vmax = na.nanmax(vals[vit])
        vals[vi] = 0.0
        if logIt:
            self.norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax, clip=False)
            self.ticker = matplotlib.ticker.LogLocator(subs=[0.1, 0.5, 1])
        else:
            self.norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=False)
            self.ticker = None
        self.cmap = matplotlib.cm.get_cmap()
        self.cmap.set_bad("k")
        self.cmap.set_under("k")
        self.cmap.set_over("k")
        self.image = self.axes.pcolor(self.x_bins, self.y_bins, \
                                      vals,shading='flat', \
                                      norm=self.norm)
        #self.ticker = matplotlib.ticker.LogLocator(subs=[0.25, 0.5, 0.75, 1])

        self.colorbar = self.figure.colorbar(self.image, \
                                             extend='neither', \
                                             shrink=0.95, cmap=self.cmap, \
                               ticks = self.ticker, format="%0.2e" )
        self.colorbar.notify(self.image)

        self.autoset_label(fields[0], self.axes.set_xlabel)
        self.autoset_label(fields[1], self.axes.set_ylabel)
        self.autoset_label(fields[2], self.colorbar.set_label)

    def generatePrefix(self, prefix):
        self.prefix = "_".join([prefix, self.typeName, \
            self.axisNames['X'], self.axisNames['Y'], \
            self.axisNames['Z']])
        self["Field1"] = self.axisNames["X"]
        self["Field2"] = self.axisNames["Y"]
        self["Field3"] = self.axisNames["Z"]

