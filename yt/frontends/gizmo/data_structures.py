import os

from yt.frontends.gadget.data_structures import GadgetHDF5Dataset
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import GizmoFieldInfo


class GizmoDataset(GadgetHDF5Dataset):
    _field_info_class = GizmoFieldInfo

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        need_groups = ["Header"]
        veto_groups = ["FOF", "Group", "Subhalo"]
        valid_fname = filename
        # If passed arg is a directory, look for the .0 file in that dir
        if os.path.isdir(filename):
            valid_files = []
            for f in os.listdir(filename):
                fname = os.path.join(filename, f)
                fext = os.path.splitext(fname)[-1]
                if (
                    (".0" in f)
                    and (fext not in {".ewah", ".kdtree"})
                    and os.path.isfile(fname)
                ):
                    valid_files.append(fname)
            if len(valid_files) != 1:
                return False
            else:
                valid_fname = valid_files[0]
        try:
            with h5py.File(valid_fname, mode="r") as fh:
                dmetal = "/PartType0/Metallicity"
                if dmetal not in fh or fh[dmetal].shape[1] < 11:
                    return False

                return all(ng in fh["/"] for ng in need_groups) and not any(
                    vg in fh["/"] for vg in veto_groups
                )

        except Exception:
            return False

    def _set_code_unit_attributes(self):
        super()._set_code_unit_attributes()
