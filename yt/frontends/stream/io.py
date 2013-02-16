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

from yt.utilities.io_handler import \
    BaseIOHandler, _axis_ids
from yt.utilities.logger import ytLogger as mylog

class IOHandlerStream(BaseIOHandler):

    _data_style = "stream"

    def __init__(self, stream_handler):
        self.fields = stream_handler.fields
        BaseIOHandler.__init__(self)

    def _read_data(self, grid, field):
        # This is where we implement processor-locking
        #if grid.id not in self.grids_in_memory:
        #    mylog.error("Was asked for %s but I have %s", grid.id, self.grids_in_memory.keys())
        #    raise KeyError
        tr = self.fields[grid.id][field].copy()
        # If it's particles, we copy.
        return tr

    def modify(self, field):
        return field

    def _read_field_names(self, grid):
        return self.fields[grid.id].keys()

    @property
    def _read_exception(self):
        return KeyError

