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
import numpy as na

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog
import cStringIO

class IOHandlerRAMSES(BaseIOHandler):
    _data_style = "ramses"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        # Chunks in this case will have affiliated domain subset objects
        # Each domain subset will contain a hydro_offset array, which gives
        # pointers to level-by-level hydro information
        n = 0
        tr = dict((f, na.empty(size, dtype='float64')) for f in fields)
        cp = 0
        for chunk in chunks:
            for subset in chunk.objs:
                # Now we read the entire thing
                f = open(subset.domain.hydro_fn, "rb")
                # This contains the boundary information, so we skim through
                # and pick off the right vectors
                content = cStringIO.StringIO(f.read())
                rv = subset.fill(content, fields)
                for ft, f in fields:
                    print "Filling %s with %s (%0.3e %0.3e) (%s:%s)" % (
                        f, subset.cell_count, rv[f].min(), rv[f].max(),
                        cp, cp+subset.cell_count)
                    tr[(ft, f)][cp:cp+subset.cell_count] = rv.pop(f)
                cp += subset.cell_count
        return tr

    def _read_data_set(self, grid, field):
        tr = na.zeros(grid.ActiveDimensions, dtype='float64')
        filled = na.zeros(grid.ActiveDimensions, dtype='int32')
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

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        return self._read_data_set(grid, field)[sl]

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
                    self.queue[g.id][field] = self._read_data_set(g, field)
                print "Clearing", field, domain
                self.ramses_tree.clear_tree(field, domain - 1)
        mylog.debug("Finished read of %s", sets)

    def modify(self, data): return data
