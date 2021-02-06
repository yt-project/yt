import os

import yt.units
from yt.frontends.gadget.data_structures import GadgetHDF5Dataset
from yt.utilities.definitions import sec_conversion
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import OWLSFieldInfo


class OWLSDataset(GadgetHDF5Dataset):
    _particle_mass_name = "Mass"
    _field_info_class = OWLSFieldInfo
    _time_readin = "Time_GYR"

    def _parse_parameter_file(self):

        # read values from header
        hvals = self._get_hvals()
        self.parameters = hvals

        # set features common to OWLS and Eagle
        self._set_owls_eagle()

        # Set time from value in header
        self.current_time = (
            hvals[self._time_readin] * sec_conversion["Gyr"] * yt.units.s
        )

    def _set_code_unit_attributes(self):
        self._set_owls_eagle_units()

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        need_groups = ["Constants", "Header", "Parameters", "Units"]
        veto_groups = [
            "SUBFIND",
            "FOF",
            "PartType0/ChemistryAbundances",
            "PartType0/ChemicalAbundances",
            "RuntimePars",
            "HashTable",
        ]
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
            valid_fname = valid_files[0]
        try:
            with h5py.File(valid_fname, mode="r") as fileh:
                for ng in need_groups:
                    if ng not in fileh["/"]:
                        return False
                for vg in veto_groups:
                    if vg in fileh["/"]:
                        return False
            return True
        except Exception:
            return False
