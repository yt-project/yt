"""
Data structures for gizmo frontend.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Britton Smith.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.data_objects.static_output import \
    ParticleFile
from yt.frontends.gadget.data_structures import \
    GadgetHDF5Dataset

from .fields import \
    GizmoFieldInfo

class GizmoDataset(GadgetHDF5Dataset):
    _file_class = ParticleFile
    _field_info_class = GizmoFieldInfo
    _particle_mass_name = "Masses"
    _suffix = ".hdf5"
