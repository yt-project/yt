#
# hippo_plot:
#   A module for interacting with HippoDraw
#
# Written by: Matthew Turk (mturk@stanford.edu) Nov 2006
# Modified:
#

from yt.raven import *

class EnzoPlot:
    def __init__(self, hierarchy, canvas, enzoHippo, offScreen):
        self.hierarchy = hierarchy
        self.canvas = canvas
        self.enzoHippo = enzoHippo
        self.offScreen = offScreen
        self.data = None
        self.im = {}
        self.im["ParameterFile"] = self.hierarchy.parameterFilename
    
    def saveImage(self, prefix, suffix):
        self.generatePrefix(prefix)
        fn = "%s.%s" % (self.prefix, suffix)
        #print "Saving to %s" % (fn)
        self.canvas.saveAsImage(self.plot, fn)
        if self.enzoHippo.httpPrefix:
            self.im["Filename"] = self.enzoHippo.httpPrefix + "/" \
                                + os.path.basename(fn)
        else:
            self.im["Filename"] = os.path.basename(fn)
        self.im["Type"] = self.typeName
        self.im["GeneratedAt"] = self.hierarchy["CurrentTimeIdentifier"]
        self.im["RunID"] = self.enzoHippo.runID
        return fn

    def getWidth(self):
        # First we iterate over the units.  This could be kind of slow.
        u=[]
        for item in self.hierarchy.units.items():
            u.append((item[1],item[0]))
        u.sort()
        for unit in u:
            print "\tWidth: %0.3e %s" % (self.width*unit[0], unit[1])

class EnzoRadialPlot(EnzoPlot):
    def makePlot(self, center, radius, unit, fields):
        time1=time.time()
        self.unit = unit
        if isinstance(unit, types.StringType):
            try:
                unit = self.hierarchy.units[unit]
            except:
                print "Unit %s not found, setting to 1.0" % (unit)
                unit = 1
        self.radius = radius/unit
        self.center = center
        self.width = radius
        #xs, ys, zs, vs = self.hierarchy.getSphere(center, radius, fields)
        self.sphere_data = lagos.EnzoSphere(self.hierarchy, self.center, self.radius, fields)
        self.data = self.sphere_data
        self.fields = fields

        self.dataLabel = []

        self.tuple = hippo.DataArray('NumArrayTuple')
        self.tuple.register()
        for i in range(len(fields)):
            field = fields[i]
            if fieldInfo.has_key(field):
                self.dataLabel.append(field + " (%s)" % (fieldInfo[field][0]))
            else:
                self.dataLabel.append(field)
            self.tuple.addColumn(field,self.data[field])
        self.field = self.fields[1]
        self.plotFromData(self.tuple)

        time2=time.time()
        print "Took %0.3e seconds for everything" % (time2-time1)

    def refresh(self):
        self.plot.setLabel('Y',self.dataLabel[1])
        if fieldInfo.has_key(self.field):
            sl = fieldInfo[self.field][2]
        else:
            sl = self.field in log_fields
        print "Setting log to %s" % (sl)
        self.plot.setLog('Y',sl)

    def plotFromData(self, dataTuple):
        # We assume we already have the data associated with the instance
        if self.data == None:
            print "Okay dumbface, it's not gonna work!  You need to have DATA!"
            print "Try calling makePlot"
            return
        self.plot =  hippo.Display(self.plotType, dataTuple, \
                    (tuple(self.fields)))
        self.plot.setLabel('X',self.dataLabel[0])
        self.plot.setLog('X',True)
        self.plot.setLabel('Y',self.dataLabel[1])
        if fieldInfo.has_key(self.field):
            sl = fieldInfo[self.field][2]
        else:
            sl = self.field in log_fields
        print "Setting log to %s" % (sl)
        self.plot.setLog('Y',sl)
        self.plot.setAspectRatio(1)

    def generatePrefix(self, prefix):
        self.prefix = prefix + "_%s" % (self.typeName)
        for field in self.fields:
            self.prefix += "_%s" % (field)
        self.im["FieldList"] = self.fields
        self.im["Fields"] = string.join(self.fields,"_")
        self.im["Unit"] = self.unit
        self.im["Width"] = self.radius * self.hierarchy[self.unit]
        print "Setting width to",self.im["Width"] * self.hierarchy[self.unit]

    def setWidth(self, width, conv):
        # In the future, we may want to consider re-generating the sphere based
        # on the new width fed here.  However, for now, that is a pretty
        # expensive operation that we will avoid.  Additionally, it would be
        # possible to get ALL of the data and then cut out the portions we
        # don't want, but that is not very memory-efficient.
        print "Not setting width of Radial Plot"
        return
        # For future reference, to regenerate the sphere we should first clear
        # the data tuple in HD's memory, erase the current data, and then
        # recall makePlot with the new radius.

    def addField(self, field):
        if self.hierarchy.conversionFactors.has_key(field):
            conv = self.hierarchy.conversionFactors[field]
        else:
            conv = 1
        if not self.tuple.has_key(field):
            self.tuple.addColumn(field,self.data[field]*conv)

    def switchField(self, field):
        self.addField(field)
        self.field = field
        if fieldInfo.has_key(field):
            self.dataLabel[1] = field + " (%s)" % (fieldInfo[field][0])
        else:
            self.dataLabel[1] = field
        rep = self.plot.getDataRep()
        rep.setAxisBinding('Y',self.field)
        self.refresh()

class EnzoRadialProfilePlot(EnzoRadialPlot):
    def __init__(self, hierarchy, canvas, enzoHippo, offScreen):
        self.typeName = "RadialProfile"
        self.numAxes = 2
        self.plotType = "Profile"
        EnzoPlot.__init__(self, hierarchy, canvas, enzoHippo, offScreen)

    def generatePrefix(self, prefix):
        self.prefix = prefix + "_%s" % (self.typeName)
        self.prefix += "_%s" % (self.field)
        self.im["FieldList"] = [self.field]
        self.im["Fields"] = self.field
        self.im["Unit"] = self.unit
        #self.im["Width"] = self.radius
        self.im["Width"] = self.radius * self.hierarchy[self.unit]
        #print "Setting width to",self.im["Width"] * self.hierarchy[self.unit]

    def makePlot(self, center, radius, unit, fields):
        self.unit = unit
        self.radius = radius
        if "Radius" in fields:
            i = fields.index("Radius")
            del fields[i]
        if "CellMass" in fields:
            i = fields.index("CellMass")
            del fields[i]
        fields = ["Radius"] + fields[:1] + ["CellMass"] + fields[1:]
        EnzoRadialPlot.makePlot(self, center, radius, unit, fields)
        self.field = self.fields[1]
        self.dataLabel = ["Radius (cm)", ""]

    def setYRange(self, minVal, maxVal):
        print "Setting range to %0.5e - %0.5e" % (minVal, maxVal)
        self.plot.setRange('Y', minVal, maxVal)

class EnzoRadialScatterPlot(EnzoRadialPlot):
    def __init__(self, hierarchy, canvas, enzoHippo, offScreen):
        self.typeName = "RadialScatter"
        self.numAxes = 2
        self.plotType = "Scatter Plot"
        EnzoPlot.__init__(self, hierarchy, canvas, enzoHippo, offScreen)

    def makePlot(self, center, radius, unit, fields):
        if "Radius" in fields:
            i = fields.index("Radius")
            del fields[i]
        fields = ["Radius"] + fields
        EnzoRadialPlot.makePlot(self, center, radius, unit, fields)

class EnzoTwoPhase(EnzoRadialPlot):
    def __init__(self, hierarchy, canvas, enzoHippo, offScreen):
        self.typeName = "TwoPhase"
        self.numAxes = 2
        self.plotType = "Color Plot"
        EnzoPlot.__init__(self, hierarchy, canvas, enzoHippo, offScreen)

    def plotFromData(self, dataTuple):
        # We assume we already have the data associated with the instance
        if self.data == None:
            print "Okay dumbface, it's not gonna work!  You need to have DATA!"
            print "Try calling makePlot"
            return
        self.plot =  hippo.Display("Color Plot", dataTuple, \
                    (tuple(self.fields)))
        for i in range(2):
            self.plot.setLabel(axis_names[i],self.dataLabel[i])
            if fieldInfo.has_key(self.fields[i]):
                sl = fieldInfo[self.fields[i]][2]
            else:
                sl = self.fields[i] in log_fields
            self.plot.setLog(axis_names[i],sl)
        self.plot.setAspectRatio(1)

class EnzoThreePhase(EnzoRadialPlot):
    def __init__(self, hierarchy, canvas, enzoHippo, offScreen):
        self.typeName = "ThreePhase"
        self.numAxes = 3
        self.plotType = "Profile 2D"
        EnzoPlot.__init__(self, hierarchy, canvas, enzoHippo, offScreen)

    def plotFromData(self, dataTuple):
        # We assume we already have the data associated with the instance
        if self.data == None:
            print "Okay dumbface, it's not gonna work!  You need to have DATA!"
            print "Try calling makePlot"
            return
        self.plot =  hippo.Display("Profile 2D", dataTuple, \
                    (self.fields[0],self.fields[1],self.fields[2]))
        for i in range(3):
            self.plot.setLabel(axis_names[i],self.dataLabel[i])
            if fieldInfo.has_key(self.fields[i]):
                sl = fieldInfo[self.fields[i]][2]
            else:
                sl = self.fields[i] in log_fields
            self.plot.setLog(axis_names[i],sl)
        self.plot.setAspectRatio(1)
        self.scaleBinWidth(0.25)

    def scaleBinWidth(self, scale):
        # We scale equally across both
        x_b = self.plot.getBinWidth('x')
        y_b = self.plot.getBinWidth('y')
        self.plot.setBinWidth('x', x_b*scale)
        self.plot.setBinWidth('y', y_b*scale)

class EnzoVM(EnzoPlot):
    def __init__(self, hierarchy, canvas, enzoHippo, offScreen):
        EnzoPlot.__init__(self, hierarchy, canvas, enzoHippo, offScreen)
        self.c = None
        self.width = 1
        self.plot = None

    def generatePrefix(self, prefix):
        self.prefix = "%s_%s_%s_%s" % (prefix, self.typeName, axis_names[self.axis], self.field)
        self.im["Axis"] = axis_names[self.axis]
        self.im["FieldList"] = [self.field]
        self.im["Field"] = self.field

    def setWidth(self, width, unit):
        self.im["Unit"] = unit
        self.im["Width"] = width
        if isinstance(unit, types.StringType):
            unit = self.hierarchy[unit]
        self.width = width / unit
        self.refreshDisplayWidth()

    def setCenter(self, c):
        self.c = c

    def refreshDisplayWidth(self, width=None):
        if width:
            self.width = width
        else:
            width = self.width
        l_edge_x = self.c[x_dict[self.axis]] - width/2.0
        r_edge_x = self.c[x_dict[self.axis]] + width/2.0
        l_edge_y = self.c[y_dict[self.axis]] - width/2.0
        r_edge_y = self.c[y_dict[self.axis]] + width/2.0
        self.plot.setRange('X', max(l_edge_x,0.0), min(r_edge_x,1.0))
        self.plot.setRange('Y', max(l_edge_y,0.0), min(r_edge_y,1.0))

    def plotFromData(self, dataTuple, cmap = "Blue-Green-Red-Yellow"):
        if self.field == "TemperatureOFF":
            cmap = ""
        # We assume we already have the data associated with the instance
        if self.data == None:
            print "Okay dumbface, it's not gonna work!  You need to have DATA!"
            print "Try calling makePlot"
            return
        print "Creating VM Display"
        if self.plot == None:
            self.plot =  hippo.Display("Variable Mesh", dataTuple, ('x','y',self.field,'dx','dx'))
        self.plot.setColorMap(cmap)
        self.plot.setAspectRatio(1)
        self.refresh()
    
    def refresh(self):
        self.plot.setAspectRatio(1)
        self.plot.setLabel('X','')
        self.plot.setLabel('Y','')
        self.plot.setAutoTicks('X',False)
        self.plot.setAutoTicks('Y',False)
        self.plot.setTicks('X',[],[])
        self.plot.setTicks('Y',[],[])
        if fieldInfo.has_key(self.field):
            sl = fieldInfo[self.field][2]
        else:
            sl = self.field in log_fields
        self.plot.setLog('Z',sl)
        self.plot.setLabel('X',axis_labels[self.axis][0])
        self.plot.setLabel('Y',axis_labels[self.axis][1])
        self.plot.setLabel('Z',self.dataLabel)

    def setZRange(self, minVal, maxVal):
        print "Setting range to %0.5e - %0.5e" % (minVal, maxVal)
        self.plot.setRange('Z', minVal, maxVal)

class EnzoVMSlice(EnzoVM):
    def makePlot(self, axis, field = "Density", center = None, slice_data = None):
        time1 = time.time()
        self.axis = axis
        self.typeName = "Slice"
        
        if (center == None) and (self.c == None):
            print "Searching for center"
            v, center = self.hierarchy.findMax('Density')
        if (center != None):
            self.c = center

        self.field = field
        
        print "Getting from field = %s at center %s on axis %s" % (field, self.c, axis)
        if slice_data == None:
            slice_data = lagos.EnzoSlice(self.hierarchy, axis, self.c[axis], field)
        
        time2 = time.time()
        print "Took %0.3e seconds to slice" % (time2-time1)
        
        self.data = slice_data
        
        self.tuple = hippo.DataArray('NumArrayTuple')
        self.tuple.register()
        
        #v1 = min(slice_data[field])
        #v2 = max(slice_data[field])
        if self.hierarchy.conversionFactors.has_key(field):
            conv = self.hierarchy.conversionFactors[field]
        else:
            conv = 1

        self.tuple.addColumn('x',self.data.x)
        self.tuple.addColumn('y',self.data.y)
        self.tuple.addColumn(field,self.data[field]*conv)
        self.tuple.addColumn('dx',self.data.dx)
        self.tuple.addColumn('dy',self.data.dy)
        
        if fieldInfo.has_key(field):
            self.dataLabel = field + " (%s)" % (fieldInfo[field][0])
        else:
            self.dataLabel = field

        #print "Min: %0.3e Max: %0.3e" % (v1, v2)
        self.plotFromData(self.tuple)
        #self.refreshDisplayWidth()
        
        time2=time.time()
        print "Took %0.3e seconds for everything" % (time2-time1)

    def addField(self, field):
        if self.hierarchy.conversionFactors.has_key(field):
            conv = self.hierarchy.conversionFactors[field]
        else:
            conv = 1
        if not self.tuple.has_key(field):
            self.tuple.addColumn(field,self.data[field]*conv)

    def switchField(self, field):
        self.addField(field)
        self.field = field
        if fieldInfo.has_key(field):
            self.dataLabel = field + " (%s)" % (fieldInfo[field][0])
        else:
            self.dataLabel = field
        rep = self.plot.getDataRep()
        rep.setAxisBinding('Z',self.field)
        self.refresh()

    def shift(self, val):
        self.data.shift(val)
        #self.tuple.dataSource().clear()
        print "Replacing x"
        self.tuple.replaceColumn('x',self.data.x)
        print "Replacing y"
        self.tuple.replaceColumn('y',self.data.y)
        print "Replacing dx"
        self.tuple.replaceColumn('dx',self.data.dx)
        print "Replacing dy"
        self.tuple.replaceColumn('dy',self.data.dy)
        for field in self.data.fields:
            print "Replacing field %s" % (field)
            self.tuple.replaceColumn(field,self.data[field])
        self.refresh()

class EnzoVMProj(EnzoVM):
    def makePlot(self, axis, field = "Density", center = None):
        time1 = time.time()
        self.axis = axis
        self.typeName = "Projection"
        
        if (center == None) and (self.c == None):
            print "Searching for center"
            v, center = self.hierarchy.findMax('Density')
        if (center != None):
            self.c = center

        self.field = field

        print "Getting from field = %s at center %s" % (field, self.c)
        projData = self.hierarchy.getProjection(axis, field)
        totalEntries = 0
        for level in projData.keys():
            totalEntries += projData[level][0].shape[0]
        x_data = zeros(totalEntries, Int64)
        y_data = zeros(totalEntries, Int64)
        z_data = zeros(totalEntries, Float64)
        dx_data = zeros(totalEntries, Float64)
        index = 0
        for level in projData.keys():
            entries = projData[level][0].shape[0]
            x_data[index:index+entries] = projData[level][0]
            y_data[index:index+entries] = projData[level][1]
            z_data[index:index+entries] = projData[level][2]
            dx_data[index:index+entries] = projData[level][3]
            index+=entries
        
        time2 = time.time()
        print "Took %0.3e seconds to project" % (time2-time1)
        
        self.data = array([(0.5+x_data)*dx_data, (0.5+y_data)*dx_data, z_data, dx_data/2.0, dx_data/2.0])
        self.data.swapaxes(0,1)

        self.tuple = hippo.DataArray('NumArrayTuple')
        self.tuple.register()
        
        v1 = self.data[:,2].min()
        v2 = self.data[:,2].max()

        self.tuple.addColumn('x',self.data[:,0])
        self.tuple.addColumn('y',self.data[:,1])
        self.tuple.addColumn(field,self.data[:,2])
        self.tuple.addColumn('dx',self.data[:,3])
        self.tuple.addColumn('dy',self.data[:,4])
        
        #for i in range(5):
            #self.tuple.addColumn(vm_axis_names[i],self.data[:,i].copy())

        if fieldInfo.has_key(field):
            self.dataLabel = field + " (%s)" % (fieldInfo[field][0])
        else:
            self.dataLabel = field

        print "Min: %0.3e Max: %0.3e" % (v1, v2)
        self.plotFromData(self.tuple)
        #self.refreshDisplayWidth()
        
        time2=time.time()
        print "Took %0.3e seconds for everything" % (time2-time1)

