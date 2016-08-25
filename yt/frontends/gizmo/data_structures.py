"""
Data structures for Gizmo frontend.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py

from yt.frontends.gadget.data_structures import \
    GadgetHDF5Dataset

from .fields import \
    GizmoFieldInfo

class GizmoDataset(GadgetHDF5Dataset):
    _field_info_class = GizmoFieldInfo

    @classmethod
    def _is_valid(self, *args, **kwargs):
        need_groups = ['Header']
        veto_groups = ['FOF', 'Group', 'Subhalo']
        valid = True
        try:
            fh = h5py.File(args[0], mode='r')
            valid = all(ng in fh["/"] for ng in need_groups) and \
              not any(vg in fh["/"] for vg in veto_groups)
            dmetal = "/PartType0/Metallicity"
            if dmetal not in fh or fh[dmetal].shape[1] not in (11, 17):
                valid = False
            fh.close()
        except:
            valid = False
            pass
        return valid
