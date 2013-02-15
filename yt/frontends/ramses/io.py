"""
RAMSES-specific IO

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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

from collections import defaultdict
import numpy as np

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

class IOHandlerRAMSES(BaseIOHandler):
    _data_style = "ramses"

    def __init__(self, ramses_tree, *args, **kwargs):
        self.ramses_tree = ramses_tree
        BaseIOHandler.__init__(self, *args, **kwargs)

    def _read_data(self, grid, field):
        tr = np.zeros(grid.ActiveDimensions, dtype='float64')
        filled = np.zeros(grid.ActiveDimensions, dtype='int32')
        to_fill = grid.ActiveDimensions.prod()
        grids = [grid]
        l_delta = 0
        varindex = self.ramses_tree.field_ind[field]
        while to_fill > 0 and len(grids) > 0:
            next_grids = []
            for g in grids:
                to_fill -= self.ramses_tree.read_grid(varindex, field,
                        grid.get_global_startindex(), grid.ActiveDimensions,
                        tr, filled, g.Level, 2**l_delta, g.locations)
                next_grids += g.Parent
            grids = next_grids
            l_delta += 1
        return tr

    def preload(self, grids, sets):
        if len(grids) == 0: return
        domain_keys = defaultdict(list)
        pf_field_list = grids[0].pf.h.field_list
        sets = [dset for dset in list(sets) if dset in pf_field_list]
        exc = self._read_exception
        for g in grids:
            domain_keys[g.domain].append(g)
        for domain, grids in domain_keys.items():
            mylog.debug("Starting read of domain %s (%s)", domain, sets)
            for field in sets:
                for g in grids:
                    self.queue[g.id][field] = self._read_data(g, field)
                print "Clearing", field, domain
                self.ramses_tree.clear_tree(field, domain - 1)
        mylog.debug("Finished read of %s", sets)

    def modify(self, data): return data
