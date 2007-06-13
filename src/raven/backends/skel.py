"""
Skeleton of a backend.  Implement all of these things, and you have a new
backend.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

engineVals = {}

def Initialize(self, **kwargs):
    pass

def CleanUp(self, **kwargs):
    pass

class RavenPlot:
    pass

class VMPlot(RavenPlot):
    pass

class ProjectionPlot(VMPlot):
    pass

class SlicePlot(VMPlot):
    pass

class LinePlot(RavenPlot):
    pass

class ProfilePlot(LinePlot):
    pass

class ACProfilePlot(LinePlot):
    pass

class RegionPlot(RavenPlot):
    pass

class TwoPhasePlot(RegionPlot):
    pass

class ThreePhasePlot(RegionPlot):
    pass
