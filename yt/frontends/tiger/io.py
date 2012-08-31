"""
TIGER-specific IO functions

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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

from yt.utilities.io_handler import \
           BaseIOHandler

class IOHandlerTiger(BaseIOHandler):
    _data_style = "tiger"
    _offset = 36

    def __init__(self, *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)
        self._memmaps = {}

    def _read_data_set(self, grid, field):
        fn = grid.pf.basename + grid.hierarchy.file_mapping[field]
        LD = np.array(grid.left_dims, dtype='int64')
        SS = np.array(grid.ActiveDimensions, dtype='int64')
        RS = np.array(grid.pf.root_size, dtype='int64')
        data = au.read_tiger_section(fn, LD, SS, RS).astype("float64")
        return data

    def _read_data_slice(self, grid, field, axis, coord):
        fn = grid.pf.basename + grid.hierarchy.file_mapping[field]
        LD = np.array(grid.left_dims, dtype='int64').copy()
        SS = np.array(grid.ActiveDimensions, dtype='int64').copy()
        RS = np.array(grid.pf.root_size, dtype='int64').copy()
        LD[axis] += coord
        SS[axis] = 1
        data = au.read_tiger_section(fn, LD, SS, RS).astype("float64")
        return data
