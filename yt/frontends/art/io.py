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
import struct
import pdb

from yt.utilities.io_handler import \
    BaseIOHandler
import numpy as na

from yt.utilities.io_handler import \
    BaseIOHandler
import yt.utilities.amr_utils as au

class IOHandlerART(BaseIOHandler):
    _data_style = "art"

    def __init__(self, filename, nhydro_vars, level_info, level_offsets,
                 *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)
        self.filename = filename
        self.nhydro_vars = nhydro_vars
        self.level_info = level_info
        self.level_offsets = level_offsets
        self.level_data = {}

    def preload_level(self, level):
        if level in self.level_data: return
        if level == 0:
            self.preload_root_level()
            return
        f = open(self.filename, 'rb')
        f.seek(self.level_offsets[level])
        ncells = 8*self.level_info[level]
        nvals = ncells * (self.nhydro_vars + 6) # 2 vars, 2 pads
        arr = na.fromfile(f, dtype='>f', count=nvals)
        arr = arr.reshape((self.nhydro_vars+6, ncells), order="F")
        arr = arr[3:-1,:].astype("float64")
        self.level_data[level] = arr

    def preload_root_level(self):
        f = open(self.filename, 'rb')
        f.seek(self.level_offsets[0] + 4) # Ditch the header
        ncells = self.level_info[0]
        #pdb.set_trace()
        nhvals = ncells * (self.nhydro_vars) # 0 vars, 0 pads
        hvar = na.fromfile(f, dtype='>f', count=nhvals).astype("float64")
        hvar = hvar.reshape((self.nhydro_vars, ncells), order="F")
        na.fromfile(f,dtype='>i',count=2) #throw away the pads
        nvars = ncells * (2) # 0 vars, 0 pads
        var = na.fromfile(f, dtype='>f', count=nvars).astype("float64")
        var = var.reshape((2, ncells), order="F")
        arr = na.concatenate((hvar,var))
        self.level_data[0] = arr

    def clear_level(self, level):
        self.level_data.pop(level, None)
        
    def _read_data_set(self, grid, field):
        pf = grid.pf
        field_id = grid.pf.h.field_list.index(field)
        if grid.Level == 0: # We only have one root grid
            self.preload_level(0)
            tr = self.level_data[0][field_id,:].reshape(
                    pf.domain_dimensions, order="F").copy()
            return tr.swapaxes(0, 2)
        tr = na.zeros(grid.ActiveDimensions, dtype='float64')
        filled = na.zeros(grid.ActiveDimensions, dtype='int32')
        to_fill = grid.ActiveDimensions.prod()
        grids = [grid]
        l_delta = 0
        while to_fill > 0 and len(grids) > 0:
            next_grids = []
            for g in grids:
                self.preload_level(g.Level)
                #print "Filling %s from %s (%s)" % (grid, g, g.Level)
                to_fill -= au.read_art_grid(field_id, 
                        grid.get_global_startindex(), grid.ActiveDimensions,
                        tr, filled, self.level_data[g.Level],
                        g.Level, 2**l_delta, g.locations)
                next_grids += g.Parent
            grids = next_grids
            l_delta += 1
        return tr

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        return self._read_data_set(grid, field)[sl]


