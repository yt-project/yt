"""
Data structures for a generic SPH/Gadget frontend.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2012 Matthew Turk.  All Rights Reserved.

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

import h5py
import numpy as np
from itertools import izip

from yt.funcs import *
from yt.geometry.oct_geometry_handler import \
    OctreeGeometryHandler
from yt.geometry.geometry_handler import \
    GeometryHandler, YTDataChunk
from yt.data_objects.static_output import \
    StaticOutput

from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc

class GadgetDomainFile(object):
    def select(self, selector):
        pass

    def count(self, selector):
        pass

class GadgetDomainSubset(object):
    def __init__(self, domain, mask, cell_count):
        self.mask = mask
        self.domain = domain
        self.oct_handler = domain.pf.h.oct_handler
        self.cell_count = cell_count
        level_counts = self.oct_handler.count_levels(
            self.domain.pf.max_level, self.domain.domain_id, mask)
        level_counts[1:] = level_counts[:-1]
        level_counts[0] = 0
        self.level_counts = np.add.accumulate(level_counts)
