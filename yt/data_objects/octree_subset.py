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

from yt.data_objects.data_containers import \
    YTFieldData, \
    YTDataContainer, \
    YTSelectionContainer
from .field_info_container import \
    NeedsGridType, \
    NeedsOriginalGrid, \
    NeedsDataField, \
    NeedsProperty, \
    NeedsParameter

class OctreeSubset(YTSelectionContainer):
    _spatial = True
    _num_ghost_zones = 0
    _num_zones = 2
    _type_name = 'octree_subset'
    _skip_add = True
    _con_args = ('domain', 'mask', 'cell_count')
    _container_fields = ("dx", "dy", "dz")

    def __init__(self, domain, mask, cell_count):
        self.field_data = YTFieldData()
        self.field_parameters = {}
        self.mask = mask
        self.domain = domain
        self.pf = domain.pf
        self.hierarchy = self.pf.hierarchy
        self.oct_handler = domain.pf.h.oct_handler
        self.cell_count = cell_count
        level_counts = self.oct_handler.count_levels(
            self.domain.pf.max_level, self.domain.domain_id, mask)
        assert(level_counts.sum() == cell_count)
        level_counts[1:] = level_counts[:-1]
        level_counts[0] = 0
        self.level_counts = np.add.accumulate(level_counts)
        self._current_particle_type = 'all'
        self._current_fluid_type = self.pf.default_fluid_type

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

    def __getitem__(self, key):
        tr = super(OctreeSubset, self).__getitem__(key)
        try:
            fields = self._determine_fields(key)
        except YTFieldTypeNotFound:
            return tr
        finfo = self.pf._get_field_info(*fields[0])
        if not finfo.particle_type:
            nz = self._num_zones + 2*self._num_ghost_zones
            n_oct = tr.shape[0] / (nz**3.0)
            dest_shape = (nz, nz, nz, n_oct)
            return tr.reshape(dest_shape)
        return tr


