"""
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


"""
Enzo's interaction with hippodraw

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from yt.raven import *
from collections import defaultdict

import hippo
import yt.lagos as lagos
from yt.funcs import *

engineVals = {}
skipAxes = ["X WIDTH", "Y WIDTH", "WEIGHT (OPTIONAL)", "DX", "DY"]

axisFieldDict = {'X':'Field1', 'Y':'Field2', 'Z':'Field3'}

def Initialize(*args, **kwargs):
    if engineVals.has_key("initialized"):
        return
    engineVals["hippoApp"] = hippo.HDApp( ) # This should throw an exception
    engineVals["canvas"] = engineVals["hippoApp"].canvas()
    engineVals["canvas"].setPlotMatrix(1,1)
    engineVals["initialized"] = True

def CleanUp(*args, **kwargs):
    engineVals["canvas"].clear()
    engineVals["canvas"].close()
    del engineVals["canvas"]
    del engineVals["hippoApp"]

class RavenPlot:
    def __init__(self, data, fields):
        self.data = data
        self.tuple = hippo.DataArray('NTuple')
        self.tuple.register()
        self.im = defaultdict(lambda: "")
        self["ParameterFile"] = self.data.hierarchy.parameterFilename
        self.axisNames = {}
        for field in fields:
            self[field] = self.data[field]
        self.hdisplay =  hippo.Display(self.plotType, self.tuple, \
                    (tuple(fields)))
        for field, axis in zip(fields, self.axes_names):
            if axis.upper() not in ["X","Y","Z"]:
                continue
            self.setupAxis(axis, field)
        self.hdisplay.setAspectRatio(1)
        engineVals["canvas"].addDisplay(self.hdisplay)

    def __setitem__(self, item, val):
        # Kind of crude, maybe not the best idea, but useful I suppose.
        # We'll need to override this for the plot types that should not be
        # converted.
        if isinstance(val, na.ndarray):
            if not self.tuple.has_key(item):
                self.tuple[item] = val * \
                    self.data.hierarchy.conversionFactors[item]
        else:
            self.im[item] = val

    def setAxisToField(self, axis, field):
        self[field] = self.data[field]
        self.hdisplay.getDataRep().setAxisBinding(axis, field)
        self.setupAxis(axis, field)

    def setupAxis(self, axis, field):
        dataLabel = field
        if lagos.fieldInfo.has_key(field):
            dataLabel += " (%s)" % (lagos.fieldInfo[field][0])
        self.hdisplay.setLabel(axis,dataLabel)
        if field in lagos.log_fields:
            self.hdisplay.setLog(axis,True)
        else:
            self.hdisplay.setLog(axis,lagos.fieldInfo[field][2])
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
        engineVals["canvas"].saveAsImage(self.hdisplay, fn)
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
        self.hdisplay.setRange('X', xmin, xmax)
        
    def set_ylim(self, ymin, ymax):
        self.hdisplay.setRange('Y', ymin, ymax)

    def set_zlim(self, zmin, zmax):
        self.hdisplay.setRange('Z', zmin, zmax)

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

class TwoPhasePlot(RavenPlot):
    def __init__(self, data, fields, width=None, unit=None):
        self.axes_names = ["X","Y","Weight (optional)"]
        self.typeName = "TwoPhase"
        self.plotType = "Color Plot"
        RavenPlot.__init__(self, data, fields)

    def generatePrefix(self, prefix):
        self.prefix = prefix + "_%s_" % (self.typeName)
        self.prefix += "_".join([self.axisNames['X'], \
                                 self.axisNames['Y']])

class ThreePhasePlot(RavenPlot):
    def __init__(self, data, fields, width=None, unit=None):
        self.axes_names = ["X","Y","Z"]
        self.typeName = "ThreePhase"
        self.plotType = "Profile 2D"
        RavenPlot.__init__(self, data, fields)

    def scaleBinWidth(self, scale):
        # We scale equally across both
        x_b = self.hdisplay.getBinWidth('X')
        y_b = self.hdisplay.getBinWidth('Y')
        self.hdisplay.setBinWidth('X', x_b*scale)
        self.hdisplay.setBinWidth('Y', y_b*scale)

    def generatePrefix(self, prefix):
        self.prefix = prefix + "_%s_" % (self.typeName)
        self.prefix += "_".join([self.axisNames['X'], \
                                 self.axisNames['Y'], \
                                 self.axisNames['Z']])

class VMPlot(RavenPlot):
    def __init__(self, data, field):
        self.axes_names = ["X","Y","Z", "X width", "Y width"]
        self.plotType = "Variable Mesh"
        fields = ['X', 'Y', field, 'X width', 'Y width']
        RavenPlot.__init__(self, data, fields)
        self.axisNames["Z"] = field
        self.checkColormap(field)
        self.set_width(1,'1')
        self.hdisplay.setLabel('X','')
        self.hdisplay.setLabel('Y','')
        self.hdisplay.setAutoTicks('X',False)
        self.hdisplay.setAutoTicks('Y',False)
        self.hdisplay.setTicks('X',[],[])
        self.hdisplay.setTicks('Y',[],[])

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
        self.hdisplay.setColorMap(cmap)

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
        self.hdisplay.setLabel('Y',"Entries per Bin")

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

