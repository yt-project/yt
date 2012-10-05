"""
Utilities for flagging zones for refinement in a dataset

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

import numpy as np # For modern purposes

flagging_method_registry = {}

class FlaggingMethod(object):
    _skip_add = False
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if hasattr(cls, "_type_name") and not cls._skip_add:
                flagging_method_registry[cls._type_name] = cls

class OverDensity(FlaggingMethod):
    _type_name = "overdensity"
    def __init__(self, over_density):
        self.over_density = over_density

    def __call__(self, grid):
        rho = grid["Density"] / (grid.pf.refine_by**grid.Level)
        return (rho > self.over_density)

class FlaggingGrid(object):
    def __init__(self, grid, methods):
        self.grid = grid
        self.sigs = []
        flagged = np.zeros(self.grid.ActiveDimensions, dtype="bool")
        for method in methods:
            flagged |= method(self.grid)
        for dim in range(3):
            d1 = (dim + 1) % 3
            d2 = (dim == 0)
            self.sigs.append(flagged.sum(axis=d1).sum(axis=d2))
        self.flagged = flagged

    def find_by_zero_signature(self):
        ge = []
        for dim in range(3):
            sig = self.sigs[dim]
            grid_ends = np.zeros((sig.size, 2))
            ng = 0
            i = 0
            while i < sig.size:
                if sig[i] != 0:
                    grid_ends[ng, 0] = i
                    while i < sig.size and sig[i] != 0:
                        i += 1
                    grid_ends[ng, 1] = i - 1
                    ng += 1
                i += 1
            ge.append(grid_ends[:ng,:])
        return ge

    def find_by_second_derivative(self):
        ze = []
        for dim in range(3):
            sig = self.sigs[dim]
            sd = sig[:-2] - 2.0*sig[1:-1] + sig[2:]
            grid_ends = np.zeros((sig.size, 2))
            ng = 0
            center = int((self.flagged.shape[dim] - 1) / 2)
            strength = zero_strength = 0
            for i in range(1, sig.size-1):
                # Note that sd is offset by one
                if sd[i-1] * sd[i] < 0:
                    strength = np.abs(sd[i-1] - sd[i])
                    if strength > zero_strength or \
                       (strength == zero_strength and np.abs(center - i) < np.abs(zero_cross -i )):
                        zero_strength = strength
                        zero_cross = i
            ze.append(zero_cross)
        return ze
