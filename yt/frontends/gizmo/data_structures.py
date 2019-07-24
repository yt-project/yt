import os

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
        valid_fname = args[0]
        # If passed arg is a directory, look for the .0 file in that dir
        if os.path.isdir(args[0]):
            valid_files = []
            for f in os.listdir(args[0]):
                fname = os.path.join(args[0], f)
                if ('.0' in f) and ('.ewah' not in f) and os.path.isfile(fname):
                    valid_files.append(fname)
            if len(valid_files) == 0:
                valid = False
            elif len(valid_files) > 1:
                valid = False
            else:
                valid_fname = valid_files[0]
        try:
            fh = h5py.File(valid_fname, mode='r')
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

    def _set_code_unit_attributes(self):
        super(GizmoDataset, self)._set_code_unit_attributes()
        self.magnetic_unit = self.quan(1.0, "gauss")
