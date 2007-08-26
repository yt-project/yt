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
All of the base-level stuff for plotting.

Think of this as a way of getting rid of EnzoHippo.  We should have access to
all of the engine-appropriate methods.

@todo: Implement regions

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from yt.raven import *

class PlotCollection:
    def __init__(self, pf, submitToDeliverator=-1):
        be.Initialize()
        self.plots = []
        self.runID = submitToDeliverator
        self.pf = pf
        if submitToDeliverator > 0:
            self.submit = True
            self.RunID = submitToDeliverator
            r=deliveration.SubmitParameterFile(\
                submitToDeliverator, self.pf)
            mylog.debug("Received response '%s'", r)
            self.httpPrefix = ytcfg["raven","httpPrefix"] % self.pf
        else:
            self.submit = False
    def save(self, basename, format="png"):
        fn = []
        for plot in self.plots:
            fn.append(plot.saveImage(basename, \
                      format="png", submit=self.runID))
            if self.submit:
                im = plot.im.copy()
                #print self.httpPrefix, "/", os.path.basename(fn[-1])
                im["Filename"] = self.httpPrefix + "/" \
                                + os.path.basename(fn[-1])
                im["RunID"] = self.RunID
                deliveration.SubmitImage(\
                      self.pf.hierarchy, im)
            mylog.info("Saved %s", fn[-1])
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
            v,center = self.pf.hierarchy.findMax("Density")
        if coord == None:
            coord = center[axis]
        slice = self.pf.hierarchy.slice(axis, coord, field, center)
        p = self.addPlot(be.SlicePlot(slice, field))
        p["Axis"] = lagos.axis_names[axis]
        return p
    def addProjection(self, field, axis, weightField=None, center=None):
        # Make the projection object here
        # Pass it in to the engine to get a SlicePlot, which we then append
        # onto self.plots
        if center == None:
            v, center = self.pf.hierarchy.findMax("Density")
        proj = self.pf.hierarchy.proj(axis, field, weightField, center=center)
        p = self.addPlot(be.ProjectionPlot(proj, field))
        p["Axis"] = lagos.axis_names[axis]
        return p
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
            v,center = self.pf.hierarchy.findMax("Density")
        r = radius/self.pf[unit]
        sphere = self.pf.hierarchy.sphere(center, r, fields)
        p = self.addPlot(be.TwoPhasePlot(sphere, fields, width=radius, unit=unit))
        p["Width"] = radius
        p["Unit"] = unit
        p["Axis"] = None
        return p 
    def addThreePhaseSphere(self, radius, unit, fields, center=None):
        if center == None:
            v,center = self.pf.hierarchy.findMax("Density")
        r = radius/self.pf[unit]
        sphere = self.pf.hierarchy.sphere(center, r, fields)
        p = self.addPlot(be.ThreePhasePlot(sphere, fields, width=radius, unit=unit))
        p["Width"] = radius
        p["Unit"] = unit
        p["Axis"] = None
        return p
