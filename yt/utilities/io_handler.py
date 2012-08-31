"""
The data-file handling functions

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

import yt.utilities.lib as au
import exceptions
import cPickle

_axis_ids = {0:2,1:1,2:0}

io_registry = {}

class BaseIOHandler(object):

    _data_style = None
    _particle_reader = False

    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if hasattr(cls, "_data_style"):
                print 'Registering Class ', cls ,' with datastyle ', cls._data_style
                io_registry[cls._data_style] = cls

    def __init__(self):
        self.queue = defaultdict(dict)

    # We need a function for reading a list of sets
    # and a function for *popping* from a queue all the appropriate sets

    def preload(self, grids, sets):
        pass

    def pop(self, grid, field):
        if grid.id in self.queue and field in self.queue[grid.id]:
            return self.modify(self.queue[grid.id].pop(field))
        else:
            # We only read the one set and do not store it if it isn't pre-loaded
            return self._read_data_set(grid, field)

    def peek(self, grid, field):
        return self.queue[grid.id].get(field, None)

    def push(self, grid, field, data):
        if grid.id in self.queue and field in self.queue[grid.id]:
            raise ValueError
        self.queue[grid][field] = data

    # Now we define our interface
    def _read_data_set(self, grid, field):
        pass

    def _read_data_slice(self, grid, field, axis, coord):
        pass

    def _read_field_names(self, grid):
        pass

    @property
    def _read_exception(self):
        return None

class IOHandlerExtracted(BaseIOHandler):

    _data_style = 'extracted'

    def _read_data_set(self, grid, field):
        return (grid.base_grid[field] / grid.base_grid.convert(field))

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        return grid.base_grid[field][tuple(sl)] / grid.base_grid.convert(field)
