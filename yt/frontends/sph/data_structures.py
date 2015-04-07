"""
Data structures for SPH frontends.




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.data_objects.static_output import \
    Dataset

class ParticleDataset(Dataset):
    _unit_base = None
    over_refine_factor = 1
    filter_bbox = False
