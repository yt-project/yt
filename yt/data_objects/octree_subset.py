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
import yt.geometry.particle_deposit as particle_deposit

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
        self._last_mask = None
        self._last_selector_id = None
        self._current_particle_type = 'all'
        self._current_fluid_type = self.pf.default_fluid_type

    def _generate_container_field(self, field):
        if self._current_chunk is None:
            self.hierarchy._identify_base_chunk(self)
        if field == "dx":
            return self._current_chunk.fwidth[:,0]
        elif field == "dy":
            return self._current_chunk.fwidth[:,1]
        elif field == "dz":
            return self._current_chunk.fwidth[:,2]

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
            # We may need to reshape the field, if it is being queried from
            # field_data.  If it's already cached, it just passes through.
            if len(tr.shape) < 4:
                tr = self._reshape_vals(tr)
            return tr
        return tr

    def _reshape_vals(self, arr):
        nz = self._num_zones + 2*self._num_ghost_zones
        n_oct = arr.shape[0] / (nz**3.0)
        arr = arr.reshape((nz, nz, nz, n_oct), order="F")
        return arr

    _domain_ind = None

    @property
    def domain_ind(self):
        if self._domain_ind is None:
            di = self.oct_handler.domain_ind(self.mask, self.domain.domain_id)
            self._domain_ind = di
        return self._domain_ind

    def deposit(self, positions, fields = None, method = None):
        # Here we perform our particle deposition.
        cls = getattr(particle_deposit, "deposit_%s" % method, None)
        if cls is None:
            raise YTParticleDepositionNotImplemented(method)
        nvals = (self.domain_ind >= 0).sum() * 8
        op = cls(nvals) # We allocate number of zones, not number of octs
        op.initialize()
        op.process_octree(self.oct_handler, self.domain_ind, positions, fields,
                          self.domain.domain_id)
        vals = op.finalize()
        return self._reshape_vals(vals)

    def select(self, selector):
        if id(selector) == self._last_selector_id:
            return self._last_mask
        self._last_mask = self.oct_handler.domain_mask(
                self.mask, self.domain.domain_id)
        if self._last_mask.sum() == 0: return None
        self._last_selector_id = id(selector)
        return self._last_mask

    def count(self, selector):
        if id(selector) == self._last_selector_id:
            if self._last_mask is None: return 0
            return self._last_mask.sum()
        self.select(selector)
        return self.count(selector)

    def count_particles(self, selector, x, y, z):
        # We don't cache the selector results
        count = selector.count_points(x,y,z)
        return count

    def select_particles(self, selector, x, y, z):
        mask = selector.select_points(x,y,z)
        return mask

class ParticleOctreeSubset(OctreeSubset):
    # Subclassing OctreeSubset is somewhat dubious.
    # This is some subset of an octree.  Note that the sum of subsets of an
    # octree may multiply include data files.  While we can attempt to mitigate
    # this, it's unavoidable for many types of data storage on disk.
    _type_name = 'particle_octree_subset'
    _con_args = ('data_files', 'pf', 'min_ind', 'max_ind')
    def __init__(self, base_selector, data_files, pf, min_ind = 0, max_ind = 0):
        # The first attempt at this will not work in parallel.
        self.data_files = data_files
        self.field_data = YTFieldData()
        self.field_parameters = {}
        self.pf = pf
        self.hierarchy = self.pf.hierarchy
        self.oct_handler = pf.h.oct_handler
        self.min_ind = min_ind
        self.max_ind = max_ind
        if max_ind == 0: max_ind = (1 << 63)
        self._last_mask = None
        self._last_selector_id = None
        self._current_particle_type = 'all'
        self._current_fluid_type = self.pf.default_fluid_type
        self.base_selector = base_selector
    
    _domain_ind = None

    @property
    def domain_ind(self):
        if self._domain_ind is None:
            mask = self.selector.select_octs(self.oct_handler)
            di = self.oct_handler.domain_ind(mask)
            self._domain_ind = di
        return self._domain_ind

    def deposit(self, positions, fields = None, method = None):
        # Here we perform our particle deposition.
        cls = getattr(particle_deposit, "deposit_%s" % method, None)
        if cls is None:
            raise YTParticleDepositionNotImplemented(method)
        nvals = (self.domain_ind >= 0).sum() * 8
        op = cls(nvals) # We allocate number of zones, not number of octs
        op.initialize()
        op.process_octree(self.oct_handler, self.domain_ind, positions, fields, 0)
        vals = op.finalize()
        return self._reshape_vals(vals)
    
    def select_icoords(self, dobj):
        return self.oct_handler.icoords(dobj.selector)

    def select_fcoords(self, dobj):
        return self.oct_handler.fcoords(dobj.selector)

    def select_fwidth(self, dobj):
        # Recall domain_dimensions is the number of cells, not octs
        base_dx = (self.pf.domain_width /
                   self.pf.domain_dimensions)
        dds = (2**self.select_ires(dobj))
        widths = np.empty((dds.shape[0], 3), dtype="float64")
        for i in range(3):
            widths[:,i] = base_dx[i] / dds
        return widths

    def select_ires(self, dobj):
        return self.oct_handler.ires(dobj.selector)

    def select(self, selector):
        if id(selector) == self._last_selector_id:
            return self._last_mask
        self._last_mask = self.oct_handler.domain_mask(
                self.selector)
        if self._last_mask.sum() == 0: return None
        self._last_selector_id = id(selector)
        return self._last_mask

    def count(self, selector):
        if id(selector) == self._last_selector_id:
            if self._last_mask is None: return 0
            return self._last_mask.sum()
        self.select(selector)
        return self.count(selector)

    def count_particles(self, selector, x, y, z):
        # We don't cache the selector results
        count = selector.count_points(x,y,z)
        return count

    def select_particles(self, selector, x, y, z):
        mask = selector.select_points(x,y,z)
        return mask

