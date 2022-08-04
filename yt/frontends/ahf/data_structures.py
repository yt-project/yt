import glob
import os

import numpy as np

from yt.data_objects.static_output import Dataset
from yt.frontends.halo_catalog.data_structures import HaloCatalogFile
from yt.funcs import setdefaultattr
from yt.geometry.particle_geometry_handler import ParticleIndex
from yt.utilities.cosmology import Cosmology

from .fields import AHFHalosFieldInfo


class AHFHalosFile(HaloCatalogFile):
    def __init__(self, ds, io, filename, file_id, range=None):
        root, _ = os.path.splitext(filename)
        candidates = glob.glob(root + "*.AHF_halos")
        if len(candidates) == 1:
            filename = candidates[0]
        else:
            raise ValueError("Too many AHF_halos files.")
        self.col_names = self._read_column_names(filename)
        super().__init__(ds, io, filename, file_id, range)

    def read_data(self, usecols=None):
        return np.genfromtxt(self.filename, names=self.col_names, usecols=usecols)

    def _read_column_names(self, filename):
        with open(filename) as f:
            line = f.readline()
            # Remove leading '#'
            line = line[1:]
            names = line.split()
            # Remove trailing '()'
            names = [name.split("(")[0] for name in names]
            return names

    def _read_particle_positions(self, ptype, f=None):
        """
        Read all particle positions in this file.
        """

        halos = self.read_data(usecols=["Xc", "Yc", "Zc"])
        pos = np.empty((halos.size, 3), dtype="float64")
        for i, ax in enumerate("XYZ"):
            pos[:, i] = halos[f"{ax}c"].astype("float64")

        return pos


class AHFHalosDataset(Dataset):
    _index_class = ParticleIndex
    _file_class = AHFHalosFile
    _field_info_class = AHFHalosFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="ahf",
        n_ref=16,
        num_zones=2,
        units_override=None,
        unit_system="cgs",
        hubble_constant=1.0,
    ):
        root, _ = os.path.splitext(filename)
        self.log_filename = root + ".log"
        self.hubble_constant = hubble_constant

        self.n_ref = n_ref
        self.num_zones = num_zones
        super().__init__(
            filename,
            dataset_type=dataset_type,
            units_override=units_override,
            unit_system=unit_system,
        )

    def _set_code_unit_attributes(self):
        setdefaultattr(self, "length_unit", self.quan(1.0, "kpccm/h"))
        setdefaultattr(self, "mass_unit", self.quan(1.0, "Msun/h"))
        setdefaultattr(self, "time_unit", self.quan(1.0, "s"))
        setdefaultattr(self, "velocity_unit", self.quan(1.0, "km/s"))

    def _parse_parameter_file(self):
        # Read all parameters.
        simu = self._read_log_simu()
        param = self._read_parameter()

        # Set up general information.
        self.filename_template = self.parameter_filename
        self.file_count = 1
        self.parameters.update(param)
        self.particle_types = "halos"
        self.particle_types_raw = "halos"

        # Set up geometrical information.
        self.refine_by = 2
        self.dimensionality = 3
        nz = self.num_zones
        self.domain_dimensions = np.ones(self.dimensionality, "int32") * nz
        self.domain_left_edge = np.array([0.0, 0.0, 0.0])
        # Note that boxsize is in Mpc but particle positions are in kpc.
        self.domain_right_edge = np.array([simu["boxsize"]] * 3) * 1000
        self._periodicity = (True, True, True)

        # Set up cosmological information.
        self.cosmological_simulation = 1
        self.current_redshift = param["z"]
        self.omega_lambda = simu["lambda0"]
        self.omega_matter = simu["omega0"]
        cosmo = Cosmology(
            hubble_constant=self.hubble_constant,
            omega_matter=self.omega_matter,
            omega_lambda=self.omega_lambda,
        )
        self.current_time = cosmo.lookback_time(param["z"], 1e6).in_units("s")

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if not filename.endswith(".parameter"):
            return False
        with open(filename) as f:
            if f.readlines()[11].startswith("AHF"):
                return True
        return False

    # Helper methods

    def _read_log_simu(self):
        simu = {}
        with open(self.log_filename) as f:
            for l in f:
                if l.startswith("simu."):
                    name, val = l.split(":")
                    key = name.strip().split(".")[1]
                    try:
                        val = float(val)
                    except Exception:
                        val = float.fromhex(val)
                    simu[key] = val
        return simu

    def _read_parameter(self):
        param = {}
        with open(self.parameter_filename) as f:
            for l in f:
                words = l.split()
                if len(words) == 2:
                    key, val = words
                    try:
                        val = float(val)
                        param[key] = val
                    except Exception:
                        pass
        return param

    @property
    def _skip_cache(self):
        return True
