"""
Skeleton of a backend.  Implement all of these things, and you have a new
backend.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
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
