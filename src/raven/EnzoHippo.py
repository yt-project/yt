#
# hippo_plot:
#   A module for interacting with HippoDraw
#
# Written by: Matthew Turk (mturk@stanford.edu) Nov 2006
# Modified:
#

from yt.raven import *
import hippo
import yt.deliverator

class EnzoHippo:
    def __init__(self, hierarchy, app=None, canvas=None, offScreen = False,\
                 gallery = None, httpPrefix=None, submitToDeliverator = -1):
        self.hierarchy = hierarchy
        self.offScreen = True

        if app == None:
            if offScreen == True:
                print "Creating non-threaded app"
                self.app = hippo.HDApp(1)
            else:
                self.app = hippo.HDApp( )
        else:
            self.app = app

        print "App: %s" % (self.app)

        if canvas == None:
            if offScreen == True:
                print "Creating independent canvas"
                self.canvas = hippo.Canvas()
            else:
                self.canvas = self.app.canvas()
        else:
            self.canvas = canvas

        self.canvas.setPlotMatrix(1,1)

        print "Canvas: %s" % (self.canvas)

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
            print "Received response '%s'" % (r)
        else:
            self.submit = False
            
        print "Returning from init"

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
        print axis
        if axis == None:
            axis = [0,1,2]
        else:
            if isinstance(axis, types.IntType):
                axis = [axis]
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        elif len(center) != 3:
            print "Center must be a 3-tuple! Using maximum density."
            v, center = self.hierarchy.findMax('Density')
        startI = len(self.plots)
        for ax in axis:
            self.plots.append(EnzoVMSlice(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
            self.plots[-1].makePlot(ax, field, center)
            self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[startI:]

    def addProj(self, field = "Density", axis = None):
        print axis
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
            print "Plot %s" % (i)
            self.plots[i].getWidth()

    def saveImages(self, prefix, suffix='png', plotIs = None):
        if plotIs:
            if isinstance(plotIs, types.IntType):
               plotIs = [plotIs]
        else:
            plotIs = range(len(self.plots))
        fn = []
        for i in plotIs:
            print "Saving..."
            fn.append(self.plots[i].saveImage(prefix, suffix))
            print fn[-1]
            if self.gallery != None:
                self.gallery.add(self.plots[i].im)
            if self.submit == True:
                yt.deliverator.SubmitImage(\
                    self.hierarchy, self.plots[i].im)
        return fn
