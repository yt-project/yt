from typing import Type

import numpy as np

import yt.units
from yt.fields.field_info_container import FieldInfoContainer
from yt.frontends.gadget.data_structures import GadgetHDF5Dataset
from yt.frontends.owls.fields import OWLSFieldInfo
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import EagleNetworkFieldInfo


class EagleDataset(GadgetHDF5Dataset):
    _particle_mass_name = "Mass"
    _field_info_class: Type[FieldInfoContainer] = OWLSFieldInfo
    _time_readin_ = "Time"

    def _parse_parameter_file(self):

        # read values from header
        hvals = self._get_hvals()
        self.parameters = hvals

        # set features common to OWLS and Eagle
        self._set_owls_eagle()

        # Set time from analytic solution for flat LCDM universe
        a = hvals["ExpansionFactor"]
        H0 = hvals["H(z)"] / hvals["E(z)"]
        a_eq = (self.omega_matter / self.omega_lambda) ** (1.0 / 3)
        t1 = 2.0 / (3.0 * np.sqrt(self.omega_lambda))
        t2 = (a / a_eq) ** (3.0 / 2)
        t3 = np.sqrt(1.0 + (a / a_eq) ** 3)
        t = t1 * np.log(t2 + t3) / H0
        self.current_time = t * yt.units.s

    def _set_code_unit_attributes(self):
        self._set_owls_eagle_units()

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        need_groups = [
            "Config",
            "Constants",
            "HashTable",
            "Header",
            "Parameters",
            "RuntimePars",
            "Units",
        ]
        veto_groups = [
            "SUBFIND",
            "PartType0/ChemistryAbundances",
            "PartType0/ChemicalAbundances",
        ]
        valid = True
        try:
            fileh = h5py.File(filename, mode="r")
            for ng in need_groups:
                if ng not in fileh["/"]:
                    valid = False
            for vg in veto_groups:
                if vg in fileh["/"]:
                    valid = False
            fileh.close()
        except Exception:
            valid = False
            pass
        return valid


class EagleNetworkDataset(EagleDataset):
    _particle_mass_name = "Mass"
    _field_info_class = EagleNetworkFieldInfo
    _time_readin = "Time"

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        try:
            fileh = h5py.File(filename, mode="r")
            if (
                "Constants" in fileh["/"].keys()
                and "Header" in fileh["/"].keys()
                and "SUBFIND" not in fileh["/"].keys()
                and (
                    "ChemistryAbundances" in fileh["PartType0"].keys()
                    or "ChemicalAbundances" in fileh["PartType0"].keys()
                )
            ):
                fileh.close()
                return True
            fileh.close()
        except Exception:
            pass
        return False
