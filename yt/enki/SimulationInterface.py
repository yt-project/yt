"""
The python interface to running a full Enzo simulation

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

def CountLinkedList(list, property):
    count = 0
    listItem = list
    while listItem != None:
        listItem = getattr(listItem, property)
        count += 1
    return count

class Simulation:
    """
    This is the simulation.
    """
    def __init__(self):
        # I think we want all of Main here.
        pass

    def run(self, nCycles = 1e30):
        # Runs the simulation through n cycles
        for i in range(nCycles):
            self.EvolveHierarchy()

    def EvolveHierarchy(self):
        # Check for stop here

        # Synchronize with Barrier
        EnzoInterface.CommunicationBarrier()
        # Calculate timestep for root grid; this sets the timestep for output,
        # as well, so compare against that.
        # EvolveLevel is a recursive call, so we only end up calling
        # that for the root grid

    def EvolveLevel(self, level):
        pass

    def _CountLevels(self):
        return CountLinkedList(self.LevelArray, 'NextLevel')



    __nLevels = None
    nLevels = property(fget=_CountLevels)