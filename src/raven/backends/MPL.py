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
import matplotlib._image
import matplotlib.colors
from matplotlib.backends.backend_agg import FigureCanvasAgg

class _AMRImage(matplotlib.image.NonUniformImage):
    _buff = None
    def make_image(self, magnification=1.0):
        if self._A is None:
            raise RuntimeError('You must first set the image array')
        if self._buff != None:
            del self._buff

        x0, y0, v_width, v_height = self.axes.viewLim.get_bounds()
        l, b, width, height = self.axes.bbox.get_bounds()
        width *= magnification
        height *= magnification
        buff = _MPL.Pixelize(self._Ax, self._Ay, 
                                    self._Adx, self._Ady,
                                    self._Adata, int(height), int(width),
                                   (x0, x0+v_width, y0, y0+v_height),
                                    )
        self.norm.autoscale(buff)
        self.norm.vmin = buff.min()
        self.norm.vmax = buff.max()
        #print "NORM", self.norm.vmin, self.norm.vmax, buff.min(), buff.max()
        if self.next_clim[0] is not None: self.norm.vmin = self.next_clim[0]
        if self.next_clim[1] is not None: self.norm.vmax = self.next_clim[1]

        Ax = (self.cmap(self.norm(buff))*255).astype(nT.UInt8)
        im = matplotlib._image.frombyte(Ax, 1)

        self.next_clim = (None, None)

        bg = matplotlib.colors.colorConverter.to_rgba(self.axes.get_frame().get_facecolor(), 0)
        im.set_bg(*bg)
        self._buff = buff
        self._A = buff
        del buff
        return im

    def set_data(self, x, y, dx, dy, d):
        x = na.asarray(x).astype(nT.Float64)
        y = na.asarray(y).astype(nT.Float64)
        dx = na.asarray(dx).astype(nT.Float64)
        dy = na.asarray(dy).astype(nT.Float64)
        d = na.asarray(d)
        ind=na.argsort(dx)
        self._A = d[ind][::-1]
        self._Adata = d[ind][::-1]
        self._Ax = x[ind][::-1]
        self._Ay = y[ind][::-1]
        self._Adx = dx[ind][::-1]
        self._Ady = dy[ind][::-1]
        self._imcache = None
        self.next_clim = (None, None)

    def set_next_clim(self, vmin, vmax):
        self.next_clim = (vmin, vmax)

    def setWidth(self, width, unit):
        pass



def amrshow(self, x, y, dx, dy, A,
           cmap = None,
           norm = None,
           aspect=None,
           interpolation=None,
           alpha=1.0,
           vmin = None,
           vmax = None,
           origin=None,
           extent=None,
           shape=None,
           filternorm=1,
           filterrad=4.0,
           imlim=None,
           **kwargs):
    """
    """

    if not self._hold: self.cla()

    if norm is not None: assert(isinstance(norm, matplotlib.colors.Normalize))
    if cmap is not None: assert(isinstance(cmap, matplotlib.colors.Colormap))
    #if aspect is None: aspect = rcParams['image.aspect']
    if aspect is None: aspect = 1.0
    self.set_aspect(aspect)
    im = _AMRImage(self, cmap, norm, extent)

    im.set_data(x, y, dx, dy, A)
    #im.set_alpha(alpha)
    self._set_artist_props(im)
    #if norm is None and shape is None:
    #    im.set_clim(vmin, vmax)
    if vmin is not None or vmax is not None:
        im.set_clim(vmin, vmax)
    else:
        im.autoscale()

    xmin, xmax, ymin, ymax = im.get_extent()

    corners = (xmin, ymin), (xmax, ymax)
    self.update_datalim(corners)
    if self._autoscaleon:
        self.set_xlim((xmin, xmax))
        self.set_ylim((ymin, ymax))
    self.images.append(im)

    return im

matplotlib.axes.Axes.amrshow = amrshow

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
    return

def CleanUp(*args, **kwargs):
    pass

class RavenPlot:
    def __init__(self, data, fields):
        self.data = data
        # With matplotlib, we don't need to copy data back and forth.
        # Unfortunately, it is also harder to plot things, as they may not have
        # a consistent interface.  Also, the plots are designed to be much more
        # *static* than in hippodraw, so switching axes is going to be a bit of
        # a pain.
        self.im = defaultdict(lambda: "")
        self["ParameterFile"] = self.data.hierarchy.parameterFilename
        self.axisNames = {}
        for field, axis in zip(fields, self.axes_names):
            if axis.upper() not in ["X","Y","Z"]:
                continue
            self.setupAxis(axis, field)

    def __setitem__(self, item, val):
        self.im[item] = val

    def setAxisToField(self, axis, field):
        self[field] = self.data[field]
        self.mdisplay.getDataRep().setAxisBinding(axis, field)
        self.setupAxis(axis, field)

    def setupAxis(self, axis, field):
        dataLabel = field
        if lagos.fieldInfo.has_key(field):
            dataLabel += " (%s)" % (lagos.fieldInfo[field][0])
        self.mdisplay.setLabel(axis,dataLabel)
        if field in lagos.log_fields:
            self.mdisplay.setLog(axis,True)
        else:
            self.mdisplay.setLog(axis,lagos.fieldInfo[field][2])
        self.axisNames[axis] = field
        self.im[axisFieldDict[axis]] = field

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
        engineVals["canvas"].saveAsImage(self.mdisplay, fn)
        self["Type"] = self.typeName
        self["GeneratedAt"] = self.data.hierarchy["CurrentTimeIdentifier"]
        # Now submit
        return fn

    def switch_x(self, field):
        self.setAxisToField('X', field)
    
    def switch_y(self, field):
        self.setAxisToField('Y', field)
    
    def switch_z(self, field):
        self.setAxisToField('Z', field)

    def switch_weight(self, field):
        self.setAxisToField('Weight (optional)', field)

    def set_xlim(self, xmin, xmax):
        self.mdisplay.setRange('X', xmin, xmax)
        
    def set_ylim(self, ymin, ymax):
        self.mdisplay.setRange('Y', ymin, ymax)

    def set_zlim(self, zmin, zmax):
        self.mdisplay.setRange('Z', zmin, zmax)

class LinePlot(RavenPlot):
    def __init__(self, data, fields):
        self.axes_names = ["X","Y"]
        self.plotType = "XY Plot"
        RavenPlot.__init__(self, data, fields)

    def generatePrefix(self, prefix):
        self.prefix = prefix + "_%s" % (self.typeName)

    def switch_z(self, *args, **kwargs):
        pass

    def set_zlim(self, *args, **kwargs):
        pass

    def addField(self, field):
        if self.hierarchy.conversionFactors.has_key(field):
            conv = self.hierarchy.conversionFactors[field]
        else:
            conv = 1
        if not self.tuple.has_key(field):
            self.tuple[field]=self.data[field]*conv

class ProfilePlot(LinePlot):
    def __init__(self, data, fields, width=None, unit=None):
        self.typeName = "RadialProfile"
        self.plotType = "Profile"
        RavenPlot.__init__(self, data, fields)

    def generatePrefix(self, prefix):
        self.prefix = prefix + "_%s" % (self.typeName)
        self.prefix += "_%s" % (self.field)

class ACProfilePlot(LinePlot):
    pass

class ScatterPlot(RavenPlot):
    def __init__(self, data, fields, width=None, unit=None):
        self.axes_names = ["X","Y"]
        self.typeName = "RadialScatter"
        self.plotType = "Scatter Plot"
        RavenPlot.__init__(self, data, fields)

class VMPlot(RavenPlot):
    def __init__(self, data, field):
        fields = ['X', 'Y', field, 'X width', 'Y width']
        RavenPlot.__init__(self, data, fields)
        self.axisNames["Z"] = field
        self.checkColormap(field)
        self.set_width(1,'1')

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
        self.mdisplay.setColorMap(cmap)

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

    def switch_y(self, *args, **kwargs):
        pass

    def switch_x(self, *args, **kwargs):
        pass

    def switch_z(self, field):
        RavenPlot.switch_z(self, field)
        self.checkColormap(field)

class SlicePlot(VMPlot):
    def __init__(self, data, fields):
        self.typeName = "Slice"
        VMPlot.__init__(self, data, fields)
        

class ProjectionPlot(VMPlot):
    def __init__(self, data, fields):
        self.typeName = "Projection"
        VMPlot.__init__(self, data, fields)

class HistogramPlot(RavenPlot):
    def __init__(self, data, field, width=0, unit=""):
        self.typeName = "Histogram"
        self.plotType = "Histogram"
        RavenPlot.__init__(self, data, [field])
        self["Width"] = float(width)
        self["Unit"] = str(unit)
        self.mdisplay.setLabel('Y',"Entries per Bin")

    def generatePrefix(self, prefix):
        self.prefix = prefix + "_%s" % (self.typeName)
        self.prefix += "_%s" % (self.field)
        self["Field1"] = self.axisNames["X"]

    def switch_y(self, *args, **kwargs):
        pass

    def set_ylim(self, *args, **kwargs):
        pass

    def switch_z(self, *args, **kwargs):
        pass

    def set_zlim(self, *args, **kwargs):
        pass

