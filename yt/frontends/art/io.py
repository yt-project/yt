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

class IOHandlerART(BaseIOHandler):
    _data_style = "art"

    def __init__(self, art_tree, *args, **kwargs):
        self.art_tree = ramses_tree
        BaseIOHandler.__init__(self, *args, **kwargs)


    def _read_data_set(self, grid, field):
        fullfieldname = 'grid_fluid_'+field
        return self.hierarchy.pf.art[fullfieldname][grid.id]

