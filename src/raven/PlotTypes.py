"""
All of the base-level stuff for plotting.

Think of this as a way of getting rid of EnzoHippo.  We should have access to
all of the engine-appropriate methods.

@todo: Implement regions

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from yt.raven import *
import yt.deliverator

class PlotCollection:
    def __init__(self, hierarchy, submitToDeliverator=-1):
        be.Initialize()
        self.plots = []
        self.runID = submitToDeliverator
        self.hierarchy = hierarchy
        if submitToDeliverator > 0:
            self.submit = True
            self.RunID = submitToDeliverator
            r=yt.deliverator.SubmitParameterFile(\
                submitToDeliverator, self.hierarchy)
            mylog.debug("Received response '%s'", r)
            self.httpPrefix = ytcfg["raven","httpPrefix"] % self.hierarchy
        else:
            self.submit = False
    def save(self, basename, format="png"):
        fn = []
        for plot in self.plots:
            fn.append(plot.saveImage(basename, \
                      format="png", submit=self.runID))
            if self.submit:
                im = plot.im.copy()
                print self.httpPrefix, "/", os.path.basename(fn[-1])
                im["Filename"] = self.httpPrefix + "/" \
                                + os.path.basename(fn[-1])
                im["RunID"] = self.RunID
                yt.deliverator.SubmitImage(\
                      self.hierarchy, im)
        return fn
    def set_xlim(self, xmin, xmax):
        for plot in self.plots:
            plot.set_xlim(xmin, xmax)
    def set_ylim(self, ymin, ymax):
        for plot in self.plots:
            plot.set_ylim(ymin, ymax)
    def set_zlim(self, zmin, zmax):
        for plot in self.plots:
            plot.set_zlim(zmin, zmax)
    def set_lim(self, lim):
        for plot in self.plots:
            plot.set_xlim(lim[:2])
            plot.set_ylim(lim[2:])
    def set_width(self, width, unit):
        for plot in self.plots:
            plot.set_width(width, unit)
    # Now we get around to adding the plots we want
    def addPlot(self, plot):
        # Accept some instance that's already been created
        # Handy for the subplot stuff that matplotlib is good at
        # And, as long as the 'plot' object is duck-typed, we should be fine
        # with it, right?
        self.plots.append(plot)
        return plot
    def addSlice(self, field, axis, coord=None, center=None):
        # Make the slice object here
        # Pass it in to the engine to get a SlicePlot, which we then append
        # onto self.plots
        if center == None:
            v,center = self.hierarchy.findMax("Density")
        if coord == None:
            coord = center[axis]
        slice = lagos.EnzoSlice(self.hierarchy, axis, coord, field, center)
        return self.addPlot(be.SlicePlot(slice, field))
    def addProjection(self, field, axis, weightField=None, center=None):
        # Make the projection object here
        # Pass it in to the engine to get a SlicePlot, which we then append
        # onto self.plots
        if center == None:
            v, center = self.hierarchy.findMax("Density")
        proj = lagos.EnzoProj(self.hierarchy, axis, field, weightField, center=center)
        return self.addPlot(be.ProjectionPlot(proj, field))
    def addACProfile(self):
        # We are given the data already
        # This is handy for plotting time-series evolution
        return None
        return self.addPlot(be.ACProfilePlot())
    def addProfile(self):
        # Make the profile here...?
        return None
        return self.addPlot(be.ProfilePlot())
    def addTwoPhaseSphere(self, radius, unit, fields, center=None):
        if center == None:
            v,center = self.hierarchy.findMax("Density")
        r = radius/self.hierarchy[unit]
        sphere = lagos.EnzoSphere(self.hierarchy, center, r, fields)
        p = self.addPlot(be.TwoPhasePlot(sphere, fields, width=radius, unit=unit))
        p["Width"] = radius
        p["Unit"] = unit
        return p 
    def addThreePhaseSphere(self, radius, unit, fields, center=None):
        if center == None:
            v,center = self.hierarchy.findMax("Density")
        r = radius/self.hierarchy[unit]
        sphere = lagos.EnzoSphere(self.hierarchy, center, r, fields)
        p = self.addPlot(be.ThreePhasePlot(sphere, fields, width=radius, unit=unit))
        p["Width"] = radius
        p["Unit"] = unit
        return p
