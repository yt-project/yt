"""
ART-specific IO

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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

import numpy as na

from yt.utilities.io_handler import \
    BaseIOHandler
import numpy as na

from yt.utilities.io_handler import \
    BaseIOHandler
import yt.utilities.amr_utils as au

class IOHandlerART(BaseIOHandler):
    _data_style = "art"

    def __init__(self, offset_root, offset_levels, nhydro_vars, ncell,
                 *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)
        self.offset_root = offset_root
        self.offset_levels = offset_levels
        self.nhydro_vars = nhydro_vars
        self.ncell= ncell
        
    def _read_data_set(self, grid, field):
        pf = grid.pf
        tr = na.zeros(grid.ActiveDimensions, dtype='float64')
        filled = na.zeros(grid.ActiveDimensions, dtype='int32')
        to_fill = grid.ActiveDimensions.prod()
        grids = [grid]
        l_delta = 0
        field_id = grid.pf.h.field_list.index(field)
        while to_fill > 0 and len(grids) > 0:
            next_grids = []
            for g in grids:
                print "Filling %s from %s" % (grid, g)
                to_fill -= au.read_art_grid(field_id, pf.ncell,
                        grid.get_global_startindex(), grid.ActiveDimensions,
                        tr, filled, g.Level, 2**l_delta, g.locations,
                        pf.parameter_filename, pf.min_level, pf.max_level,
                        self.nhydro_vars, self.offset_root, self.offset_levels,
                        pf.level_offsets)
                next_grids += g.Parent
            grids = next_grids
            l_delta += 1
        return tr

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        return self._read_data_set(grid, field)[sl]


