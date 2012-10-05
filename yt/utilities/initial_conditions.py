"""
Painting zones in a grid

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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

import numpy as np

class DataModifier(object):
    def apply(self, pf):
        for g in pf.h.grids: self(g)

class TopHatSphere(DataModifier):
    def __init__(self, radius, center, fields):
        self.radius = radius
        self.center = center
        self.fields = fields
        
    def __call__(self, grid, sub_select = None):
        r = np.zeros(grid.ActiveDimensions, dtype="float64")
        for i, ax in enumerate("xyz"):
            np.add(r, (grid[ax] - self.center[i])**2.0, r)
        np.sqrt(r, r)
        ind = (r <= self.radius)
        if sub_select is not None:
            ind &= sub_select
        for field, val in self.fields.iteritems():
            grid[field][r < self.radius] = val

class RandomFluctuation(DataModifier):
    def __init__(self, fields):
        self.fields = fields

    def __call__(self, grid, sub_select = None):
        if sub_select is None:
            sub_select = Ellipsis
        for field, mag in self.fields.iteritems():
            vals = grid[field][sub_select]
            rc = 1.0 + (np.random.random(vals.shape) - 0.5) * mag
            grid[field][sub_select] *= rc
