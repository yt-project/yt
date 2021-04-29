import glob
import os

import numpy as np

from yt.data_objects.static_output import ParticleDataset
from yt.frontends.halo_catalog.data_structures import HaloCatalogFile
from yt.funcs import setdefaultattr
from yt.geometry.particle_geometry_handler import ParticleIndex
from yt.utilities import fortran_utils as fpu
from yt.utilities.cosmology import Cosmology

from .definitions import header_dt
from .fields import RockstarFieldInfo


class RockstarBinaryFile(HaloCatalogFile):
    def __init__(self, ds, io, filename, file_id, range):
        with open(filename, "rb") as f:
            self.header = fpu.read_cattrs(f, header_dt, "=")
            self._position_offset = f.tell()
            f.seek(0, os.SEEK_END)
            self._file_size = f.tell()

        super().__init__(ds, io, filename, file_id, range)

    def _read_particle_positions(self, ptype, f=None):
        """
        Read all particle positions in this file.
        """

        if f is None:
            close = True
            f = open(self.filename, "rb")
        else:
            close = False

        pcount = self.header["num_halos"]
        pos = np.empty((pcount, 3), dtype="float64")
        f.seek(self._position_offset, os.SEEK_SET)
        halos = np.fromfile(f, dtype=self.io._halo_dt, count=pcount)
        for i, ax in enumerate("xyz"):
            pos[:, i] = halos[f"particle_position_{ax}"].astype("float64")

        if close:
            f.close()

        return pos


class RockstarDataset(ParticleDataset):
    _index_class = ParticleIndex
    _file_class = RockstarBinaryFile
    _field_info_class = RockstarFieldInfo
    _suffix = ".bin"

    def __init__(
        self,
        filename,
        dataset_type="rockstar_binary",
        units_override=None,
        unit_system="cgs",
        index_order=None,
        index_filename=None,
    ):
        super().__init__(
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
        )

    def _parse_parameter_file(self):
        with open(self.parameter_filename, "rb") as f:
            hvals = fpu.read_cattrs(f, header_dt)
            hvals.pop("unused")
        self.dimensionality = 3
        self.refine_by = 2
        prefix = ".".join(self.parameter_filename.rsplit(".", 2)[:-2])
        self.filename_template = f"{prefix}.%(num)s{self._suffix}"
        self.file_count = len(glob.glob(prefix + ".*" + self._suffix))

        # Now we can set up things we already know.
        self.cosmological_simulation = 1
        self.current_redshift = (1.0 / hvals["scale"]) - 1.0
        self.hubble_constant = hvals["h0"]
        self.omega_lambda = hvals["Ol"]
        self.omega_matter = hvals["Om"]
        cosmo = Cosmology(
            hubble_constant=self.hubble_constant,
            omega_matter=self.omega_matter,
            omega_lambda=self.omega_lambda,
        )
        self.current_time = cosmo.lookback_time(self.current_redshift, 1e6).in_units(
            "s"
        )
        self._periodicity = (True, True, True)
        self.particle_types = "halos"
        self.particle_types_raw = "halos"

        self.domain_left_edge = np.array([0.0, 0.0, 0.0])
        self.domain_right_edge = np.array([hvals["box_size"]] * 3)

        self.domain_dimensions = np.ones(3, "int32")
        self.parameters.update(hvals)

    def _set_code_unit_attributes(self):
        z = self.current_redshift
        setdefaultattr(self, "length_unit", self.quan(1.0 / (1.0 + z), "Mpc / h"))
        setdefaultattr(self, "mass_unit", self.quan(1.0, "Msun / h"))
        setdefaultattr(self, "velocity_unit", self.quan(1.0, "km / s"))
        setdefaultattr(self, "time_unit", self.length_unit / self.velocity_unit)

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if not filename.endswith(".bin"):
            return False
        with open(filename, mode="rb") as f:
            header = fpu.read_cattrs(f, header_dt)
            if header["magic"] == 18077126535843729616:
                return True
        return False
