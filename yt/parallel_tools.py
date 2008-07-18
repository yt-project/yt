"""
Parallel data mapping techniques for yt

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

import itertools
from yt.arraytypes import *

try:
    from mpi4py import MPI
    parallel_capable = True
except ImportError:
    parallel_capable = False

class GridIterator(object):
    def __init__(self, pobj):
        self.pobj = pobj
        if hasattr(pobj, '_grids') and pobj._grids is not None:
            print "I Killed"
            self._grids = pobj._grids
        else:
            print "Yo"
            self._grids = pobj._data_source._grids
        print pobj._data_source._grids, self._grids
        self.ng = len(self._grids)

    def __iter__(self):
        self.pos = 0
        return self

    def next(self):
        # We do this manually in case
        # something else asks for us.pos
        if self.pos < len(self._grids):
            self.pos += 1
            return self._grids[self.pos - 1]
        raise StopIteration

class ParallelGridIterator(GridIterator):
    """
    This takes an object, pobj, that implements ParallelAnalysisInterface,
    and then does its thing.
    """
    def __init__(self, pobj):
        GridIterator.__init__(self, pobj)
        self._offset = MPI.COMM_WORLD.rank
        self._skip = MPI.COMM_WORLD.size
        # Note that we're doing this in advance, and with a simple means
        # of choosing them; more advanced methods will be explored later.
        self.my_grid_ids = na.mgrid[self._offset:self.ng:self._skip]
        
    def __iter__(self):
        self.pobj._initialize_parallel()
        self.pos = 0
        return self

    def next(self):
        if self.pos < len(self.my_grid_ids):
            gid = self.my_grids_ids[self.pos]
            self.pos += 1
            return self._grids[gid]
        self.pobj._finalize_parallel()
        raise StopIteration

class ParallelAnalysisInterface(object):
    _grids = None

    def _get_grids(self):
        if parallel_capable and \
           ytcfg.get_boolean("yt","parallel"):
            return ParallelGridIterator(self)
        return GridIterator(self)

    def _initialize_parallel(self):
        pass

    def _finalize_parallel(self):
        pass
