"""
This is a particle container that provides no indexing information.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.data_objects.data_containers import \
    YTFieldData, \
    YTDataContainer, \
    YTSelectionContainer
from yt.funcs import *
from yt.utilities.exceptions import \
    YTNonIndexedDataContainer
from yt.data_objects.octree_subset import \
     ParticleOctreeSubset, _use_global_octree

def _non_indexed(name):
    def _func_non_indexed(self, *args, **kwargs):
        raise YTNonIndexedDataContainer(self)
    return _func_non_indexed

class ParticleContainer(YTSelectionContainer):
    _spatial = False
    _type_name = 'particle_container'
    _skip_add = True
    _con_args = ('base_region', 'data_files', 'overlap_files')

    def __init__(self, base_region, data_files, overlap_files = [], 
                 domain_id = -1):
        self.field_data = YTFieldData()
        self.field_parameters = {}
        self.data_files = ensure_list(data_files)
        self.overlap_files = ensure_list(overlap_files)
        self.ds = self.data_files[0].ds
        self._last_mask = None
        self._last_selector_id = None
        self._current_particle_type = 'all'
        # self._current_fluid_type = self.ds.default_fluid_type
        if hasattr(base_region, "base_selector"):
            self.base_selector = base_region.base_selector
            self.base_region = base_region.base_region
        else:
            self.base_region = base_region
            self.base_selector = base_region.selector
        self.domain_id = -1
            
    @property
    def selector(self):
        raise YTDataSelectorNotImplemented(self.oc_type_name)

    def select_particles(self, selector, x, y, z):
        mask = selector.select_points(x,y,z)
        return mask

    @contextlib.contextmanager
    def _expand_data_files(self):
        old_data_files = self.data_files
        old_overlap_files = self.overlap_files
        self.data_files = list(set(self.data_files + self.overlap_files))
        self.data_files.sort()
        self.overlap_files = []
        yield self
        self.data_files = old_data_files
        self.overlap_files = old_overlap_files

    select_blocks = _non_indexed('select_blocks')
    deposit = _non_indexed('deposit')
    smooth = _non_indexed('smooth')
    select_icoords = _non_indexed('select_icoords')
    select_fcoords = _non_indexed('select_fcoords')
    select_fwidth = _non_indexed('select_fwidth')
    select_ires = _non_indexed('select_ires')
    select = _non_indexed('select')
    count = _non_indexed('count')
    count_particles = _non_indexed('count_particles')
