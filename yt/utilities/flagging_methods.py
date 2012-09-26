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

def flag_cells(grid, methods):
    flagged = np.zeros(grid.ActiveDimensions, dtype="bool")
    for method in methods:
        flagged |= method(grid)
    return flagged

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

    def __call__(self, pf, grid):
        rho = grid["Density"] / (pf.refine_by**grid.Level)
        return (rho > self.over_density)
