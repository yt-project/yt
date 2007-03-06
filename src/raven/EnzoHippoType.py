"""
Enzo's interaction with hippodraw

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from EnzoPlotTypes import *
from EnzoGallery import *
from yt.raven import *
import hippo
try:
    import yt.deliverator
except:
    mylog.warning("Deliverator import failed; all deliverator actions will fail!")

class EnzoHippo:
    def __init__(self, hierarchy, app=None, canvas=None, offScreen = False,\
                 gallery = None, httpPrefix=None, submitToDeliverator = -1):
        """
        This is an instance of  a hippodraw controller designed for use with
        Enzo data.

        @param hierarchy: The particular datadump to associate it with
        @type hierarchy: L{EnzoHierarchy<lagos.EnzoHierarchyType.EnzoHierarchy>}
        @keyword app: The hippo.HDApp() to associate it with.  (Usually only used
            inside L{yt.fido}.)
        @type app: HDApp
        @keyword canvas: the hippodraw canvas we want to associate with.
            (Usually not used.)
        @type canvas: hippo.Canvas
        @keyword offScreen: doesn't work.  But hopefully, will!
        @type offScreen: boolean
        @keyword gallery: The pickled L{EnzoGallery.RavenImageDB} filename to use
        @type gallery: string
        @keyword httpPrefix: the prefix to prepend before the image filenames
            when submitting to gallery or Deliverator
        @type httpPrefix: string
        @keyword submitToDeliverator: the RunID to submit to on The Deliverator
        @type submitToDeliverator: integer
        """
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
        """
        Add a new plot of this type

        @param fields: the fields to add to the plot
        @type fields: list of strings
        @param radius: the radius to include
        @type radius: float
        @param unit: the unit the radius is listed in
        @type unit: string
        @keyword center: the center (defaults to maximum Density)
        @type center: tuple of floats
        @return: L{EnzoRadialProfilePlot<EnzoPlotTypes.EnzoRadialProfilePlot>}
        """
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        self.plots.append(EnzoRadialProfilePlot(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
        self.plots[-1].makePlot(center, radius, unit, fields)
        self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[-1]

    def addRadialScatterPlot(self, fields, radius, unit, center=None):
        """
        Add a new plot of this type

        @param fields: the fields to add to the plot
        @type fields: list of strings
        @param radius: the radius to include
        @type radius: float
        @param unit: the unit the radius is listed in
        @type unit: string
        @keyword center: the center (defaults to maximum Density)
        @type center: tuple of floats
        @return: L{EnzoRadialProfilePlot<EnzoPlotTypes.EnzoRadialProfilePlot>}
        """
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        self.plots.append(EnzoRadialScatterPlot(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
        self.plots[-1].makePlot(center, radius, unit, fields)
        self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[-1]

    def addThreePhase(self, fields, radius, unit, center=None):
        """
        Add a new plot of this type

        @param fields: the fields to add to the plot
        @type fields: list of strings
        @param radius: the radius to include
        @type radius: float
        @param unit: the unit the radius is listed in
        @type unit: string
        @keyword center: the center (defaults to maximum Density)
        @type center: tuple of floats
        @return: L{EnzoThreePhase<EnzoPlotTypes.EnzoThreePhase>}
        """
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        self.plots.append(EnzoThreePhase(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
        self.plots[-1].makePlot(center, radius, unit, fields)
        self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[-1]

    def addTwoPhase(self, fields, radius, unit, center=None):
        """
        Add a new plot of this type

        @param fields: the fields to add to the plot
        @type fields: list of strings
        @param radius: the radius to include
        @type radius: float
        @param unit: the unit the radius is listed in
        @type unit: string
        @keyword center: the center (defaults to maximum Density)
        @type center: tuple of floats
        @return: L{EnzoTwoPhase<EnzoPlotTypes.EnzoTwoPhase>}
        """
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        self.plots.append(EnzoTwoPhase(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
        self.plots[-1].makePlot(center, radius, unit, fields)
        self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[-1]

    def addSlice(self, field = "Density", axis = None, center = None, cmap = None):
        """
        Add a new plot of this type

        @param field: the field to add to the slice
        @type field: string
        @keyword axis: the axis to slice along (defaults to all three)
        @type axis: integer, list of integers
        @keyword center: the center (defaults to maximum Density)
        @type center: tuple of floats
        @keyword cmap: the colormap
        @type cmap: string
        @return: L{EnzoVMSlice<EnzoPlotTypes.EnzoVMSlice>}
        """
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
            self.plots[-1].makePlot(ax, field, center, cmap = cmap)
            self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[startI:]

    def addProj(self, field="Density", axis=None, weight=None, center=None, \
                minLevel=0, maxLevel=None):
        if axis == None:
            axis = [0,1,2]
        else:
            if isinstance(axis, types.IntType):
                axis = [axis]
        if center == None:
            v, center = self.hierarchy.findMax('Density')
        startI = len(self.plots)
        for ax in axis:
            self.plots.append(EnzoVMProj(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
            self.plots[-1].makePlot(ax, field, center, minLevel, maxLevel, weight)
            self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[startI:]

    def setCenter(self, center, plotIs = None):
        """
        Change the center of the plots

        @param center: the new width
        @keyword plotIs: the plots to change (defaults to all)
        @type plotIs: list of ints
        """
        if plotIs:
            if isinstance(plotIs, types.IntType):
               plotIs = [plotIs]
        else:
            plotIs = range(len(self.plots))
        for i in plotIs:
            self.plots[i].setCenter(center)

    def setWidth(self, width, unit, plotIs = None):
        """
        Change the width of the plots

        @param width: the new width
        @type width: float
        @param unit: the unit the width is expressed in
        @type unit: string
        @keyword plotIs: the plots to change (defaults to all)
        @type plotIs: list of ints
        """
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
        """
        Returns the widths of the plots

        @keyword plotIs: the plots to get the widths of (defaults to all)
        @type plotIs: list of ints
        """
        if plotIs:
            if isinstance(plotIs, types.IntType):
               plotIs = [plotIs]
        else:
            plotIs = range(len(self.plots))
        for i in plotIs:
            self.plots[i].getWidth()

    def saveImages(self, prefix, suffix='png', plotIs = None):
        """
        Save images from a set of plots, or all plots

        @param prefix: the prefix to prepend when creating the filenames
        @type prefix: string
        @keyword suffix: the filetype
        @type suffix: string
        @keyword plotIs: the plots to save
        @type plotIs: list of ints
        @return: list of strings
        """
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
