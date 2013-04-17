"""
Subsets of octrees

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Matthew Turk.  All Rights Reserved.

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

class OctreeSubset(object):
    def __init__(self, domain, mask, cell_count):
        self.mask = mask
        self.domain = domain
        self.oct_handler = domain.pf.h.oct_handler
        self.cell_count = cell_count
        level_counts = self.oct_handler.count_levels(
            self.domain.pf.max_level, self.domain.domain_id, mask)
        assert(level_counts.sum() == cell_count)
        level_counts[1:] = level_counts[:-1]
        level_counts[0] = 0
        self.level_counts = np.add.accumulate(level_counts)

    def select_icoords(self, dobj):
        return self.oct_handler.icoords(self.domain.domain_id, self.mask,
                                        self.cell_count,
                                        self.level_counts.copy())

    def select_fcoords(self, dobj):
        return self.oct_handler.fcoords(self.domain.domain_id, self.mask,
                                        self.cell_count,
                                        self.level_counts.copy())

    def select_fwidth(self, dobj):
        # Recall domain_dimensions is the number of cells, not octs
        base_dx = (self.domain.pf.domain_width /
                   self.domain.pf.domain_dimensions)
        widths = np.empty((self.cell_count, 3), dtype="float64")
        dds = (2**self.select_ires(dobj))
        for i in range(3):
            widths[:,i] = base_dx[i] / dds
        return widths

    def select_ires(self, dobj):
        return self.oct_handler.ires(self.domain.domain_id, self.mask,
                                     self.cell_count,
                                     self.level_counts.copy())

