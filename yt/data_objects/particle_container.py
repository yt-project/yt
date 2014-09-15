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

def _non_indexed(self, *args, **kwargs):
    raise YTNonIndexedDataContainer(self)

class ParticleContainer(YTSelectionContainer):
    _spatial = False
    _type_name = 'particle_container'
    _skip_add = True
    _con_args = ('base_region', 'data_files')

    def __init__(self, base_region, data_files):
        self.field_data = YTFieldData()
        self.data_files = ensure_list(data_files)
        self.field_parameters = {}
        self.ds = self.data_files[0].ds
        self._current_particle_type = 'all'
        self.base_region = base_region
        self.base_selector = base_region.selector

    def select_particles(self, selector, x, y, z):
        mask = selector.select_points(x,y,z)
        return mask

    select_blocks = _non_indexed
    deposit = _non_indexed
    smooth = _non_indexed
    select_icoords = _non_indexed
    select_fcoords = _non_indexed
    select_fwidth = _non_indexed
    select_ires = _non_indexed
    select = _non_indexed
    count = _non_indexed
    count_particles = _non_indexed
