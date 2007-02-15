#
# hippo_plot:
#   A module for interacting with HippoDraw
#
# Written by: Matthew Turk (mturk@stanford.edu) Nov 2006
# Modified:
#

from yt.raven import *
import hippo
try:
    import yt.deliverator
except:
    mylog.warning("Deliverator import failed; all deliverator actions will fail!")

class EnzoHippo:
    def __init__(self, hierarchy, app=None, canvas=None, offScreen = False,\
                 gallery = None, httpPrefix=None, submitToDeliverator = -1):
        self.hierarchy = hierarchy
        self.offScreen = True

        if app == None:
            if offScreen == True:
                mylog.info("Creating non-threaded app")
                self.app = hippo.HDApp(1)
            else:
                self.app = hippo.HDApp( )
        else:
            self.app = app

        if canvas == None:
            if offScreen == True:
                mylog.info("Creating independent canvas")
                self.canvas = hippo.Canvas()
            else:
                self.canvas = self.app.canvas()
        else:
            self.canvas = canvas

        self.canvas.setPlotMatrix(1,1)

        self.httpPrefix = httpPrefix

        self.plots = []
        if gallery:
            self.gallery = RavenImageDB(gallery)
        else:
            self.gallery = None

        self.runID = submitToDeliverator
        if submitToDeliverator >= 0:
            self.submit = True
            r=yt.deliverator.SubmitParameterFile(\
                submitToDeliverator, self.hierarchy)
            mylog.info("Received response '%s'", r)
        else:
            self.submit = False
            
    def __del__(self):
        """
        Note that we delete the canvas, but not the app.
        """
        self.canvas.clear()
        self.canvas.close()

    def addRadialProfilePlot(self, fields, radius, unit, center=None):
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        self.plots.append(EnzoRadialProfilePlot(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
        self.plots[-1].makePlot(center, radius, unit, fields)
        self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[-1]

    def addRadialScatterPlot(self, fields, radius, unit, center=None):
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        self.plots.append(EnzoRadialScatterPlot(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
        self.plots[-1].makePlot(center, radius, unit, fields)
        self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[-1]

    def addThreePhase(self, fields, radius, unit, center=None):
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        self.plots.append(EnzoThreePhase(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
        self.plots[-1].makePlot(center, radius, unit, fields)
        self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[-1]

    def addTwoPhase(self, fields, radius, unit, center=None):
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        self.plots.append(EnzoTwoPhase(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
        self.plots[-1].makePlot(center, radius, unit, fields)
        self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[-1]

    def addSlice(self, field = "Density", axis = None, center = None):
        if axis == None:
            axis = [0,1,2]
        else:
            if isinstance(axis, types.IntType):
                axis = [axis]
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        elif len(center) != 3:
            mylog.warning("Center must be a 3-tuple! Using maximum density.")
            v, center = self.hierarchy.findMax('Density')
        startI = len(self.plots)
        for ax in axis:
            self.plots.append(EnzoVMSlice(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
            self.plots[-1].makePlot(ax, field, center)
            self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[startI:]

    def addProj(self, field = "Density", axis = None):
        if axis == None:
            axis = [0,1,2]
        else:
            if isinstance(axis, types.IntType):
                axis = [axis]
        v, center = self.hierarchy.findMax('Density')
        startI = len(self.plots)
        for ax in axis:
            self.plots.append(EnzoVMProj(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
            self.plots[-1].makePlot(ax, field, center)
            self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[startI:]

    def setCenter(self, center, plotIs = None):
        if plotIs:
            if isinstance(plotIs, types.IntType):
               plotIs = [plotIs]
        else:
            plotIs = range(len(self.plots))
        for i in plotIs:
            self.plots[i].setCenter(center)

    def setWidth(self, width, unit, plotIs = None):
        if plotIs:
            if isinstance(plotIs, types.IntType):
               plotIs = [plotIs]
        else:
            plotIs = range(len(self.plots))
        #if isinstance(unit, types.StringType):
            #try:
                #unit = self.hierarchy.units[unit]
            #except:
                #print "Unit %s not found, setting to 1.0" % (unit)
                #unit = 1
        for i in plotIs:
            self.plots[i].setWidth(width, unit)

    def getWidth(self, plotIs = None):
        if plotIs:
            if isinstance(plotIs, types.IntType):
               plotIs = [plotIs]
        else:
            plotIs = range(len(self.plots))
        for i in plotIs:
            self.plots[i].getWidth()

    def saveImages(self, prefix, suffix='png', plotIs = None):
        if plotIs:
            if isinstance(plotIs, types.IntType):
               plotIs = [plotIs]
        else:
            plotIs = range(len(self.plots))
        fn = []
        for i in plotIs:
            fn.append(self.plots[i].saveImage(prefix, suffix))
            mylog.info("Saved %s", fn[-1])
            if self.gallery != None:
                self.gallery.add(self.plots[i].im)
            if self.submit == True:
                yt.deliverator.SubmitImage(\
                    self.hierarchy, self.plots[i].im)
        return fn
