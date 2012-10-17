"""
Enzo-specific IO functions

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

import exceptions
import os
import numpy as np

from yt.utilities.io_handler import \
    BaseIOHandler, _axis_ids
from yt.utilities.logger import ytLogger as mylog

class IOHandlerStream(BaseIOHandler):

    _data_style = "stream"

    def __init__(self, stream_handler):
        self.fields = stream_handler.fields
        BaseIOHandler.__init__(self)

    def _read_data_set(self, grid, field):
        # This is where we implement processor-locking
        #if grid.id not in self.grids_in_memory:
        #    mylog.error("Was asked for %s but I have %s", grid.id, self.grids_in_memory.keys())
        #    raise KeyError
        tr = self.fields[grid.id][field]
        # If it's particles, we copy.
        if len(tr.shape) == 1: return tr.copy()
        # New in-place unit conversion breaks if we don't copy first
        return tr

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype != "gas" for ftype, fname in fields)):
            raise NotImplementedError
        rv = {}
        for field in fields:
            ftype, fname = field
            rv[field] = np.empty(size, dtype="float64")
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [f2 for f1, f2 in fields], ng)
        for field in fields:
            ftype, fname = field
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    mask = g.select(selector) # caches
                    if mask is None: continue
                    ds = self.fields[g.id][fname]
                    data = ds[mask]
                    rv[field][ind:ind+data.size] = data
                    ind += data.size
        return rv

    @property
    def _read_exception(self):
        return KeyError

