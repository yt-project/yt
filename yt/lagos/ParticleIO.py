"""
The particle-IO handler

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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

from yt.funcs import *
from yt.lagos import *

class ParticleIOHandler(AMR3DData):
    def __init__(self, pf, source):
        self.pf = pf
        self.data = {}
        self.source = source

    def __getitem__(self, key):
        if key not in self.data:
            self.get_data(key)
        return self.data[key]

    def get_data(self, fields):
        fields = ensure_list(fields)
        self.pf.h.io._read_particles(self.source._grids,
                fields, self._dispatch)

class ParticleIOHandlerRegion(ParticleIOHandler):
    def __init__(self, pf, source, left_edge, right_edge, periodic = False):
        self.left_edge = left_edge
        self.right_edge = right_edge
        self.periodic = periodic
        ParticleIOHandler.__init__(self, pf, source)

    def dispatch(self, grid):
        """
        The purpose of this function is to determine whether or not
        we should read the entire grid.
        """
        if self.source._is_fully_enclosed(grid):
            return True
        return False
