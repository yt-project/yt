import glob
import os
from functools import cached_property
from typing import Any, Optional

import numpy as np

from yt.data_objects.static_output import ParticleDataset
from yt.frontends.halo_catalog.data_structures import HaloCatalogFile
from yt.funcs import setdefaultattr
from yt.geometry.particle_geometry_handler import ParticleIndex
from yt.utilities import fortran_utils as fpu
from yt.utilities.cosmology import Cosmology
from yt.utilities.exceptions import YTFieldNotFound

from .definitions import header_dt
from .fields import RockstarFieldInfo


class RockstarBinaryFile(HaloCatalogFile):
    header: dict
    _position_offset: int
    _member_offset: int
    _Npart: "np.ndarray[Any, np.dtype[np.int64]]"
    _ids_halos: list[int]
    _file_size: int

    def __init__(self, ds, io, filename, file_id, range):
        with open(filename, "rb") as f:
            self.header = fpu.read_cattrs(f, header_dt, "=")
            self._position_offset = f.tell()
            pcount = self.header["num_halos"]

            halos = np.fromfile(f, dtype=io._halo_dt, count=pcount)
            self._member_offset = f.tell()
            self._ids_halos = list(halos["particle_identifier"])
            self._Npart = halos["num_p"]

            f.seek(0, os.SEEK_END)
            self._file_size = f.tell()

        expected_end = self._member_offset + 8 * self._Npart.sum()
        if expected_end != self._file_size:
            raise RuntimeError(
                f"File size {self._file_size} does not match expected size {expected_end}."
            )

        super().__init__(ds, io, filename, file_id, range)

    def _read_member(
        self, ihalo: int
    ) -> Optional["np.ndarray[Any, np.dtype[np.int64]]"]:
        if ihalo not in self._ids_halos:
            return None

        ind_halo = self._ids_halos.index(ihalo)

        ipos = self._member_offset + 8 * self._Npart[:ind_halo].sum()

        with open(self.filename, "rb") as f:
            f.seek(ipos, os.SEEK_SET)
            ids = np.fromfile(f, dtype=np.int64, count=self._Npart[ind_halo])
            return ids

    def _read_particle_positions(self, ptype: str, f=None):
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


class RockstarIndex(ParticleIndex):
    def get_member(self, ihalo: int):
        for df in self.data_files:
            members = df._read_member(ihalo)
            if members is not None:
                return members

        raise RuntimeError(f"Could not find halo {ihalo} in any data file.")


class RockstarDataset(ParticleDataset):
    _index_class = RockstarIndex
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
    def _is_valid(cls, filename: str, *args, **kwargs) -> bool:
        if not filename.endswith(".bin"):
            return False
        try:
            with open(filename, mode="rb") as f:
                header = fpu.read_cattrs(f, header_dt)
        except OSError:
            return False
        else:
            return header["magic"] == 18077126535843729616

    def halo(self, ptype, particle_identifier):
        return RockstarHaloContainer(
            ptype,
            particle_identifier,
            parent_ds=None,
            halo_ds=self,
        )


class RockstarHaloContainer:
    def __init__(self, ptype, particle_identifier, *, parent_ds, halo_ds):
        if ptype not in halo_ds.particle_types_raw:
            raise RuntimeError(
                f'Possible halo types are {halo_ds.particle_types_raw}, supplied "{ptype}".'
            )

        self.ds = parent_ds
        self.halo_ds = halo_ds
        self.ptype = ptype
        self.particle_identifier = particle_identifier

    def __repr__(self):
        return f"{self.halo_ds}_{self.ptype}_{self.particle_identifier:09d}"

    def __getitem__(self, key):
        if isinstance(key, tuple):
            ptype, field = key
        else:
            ptype = self.ptype
            field = key

        data = {
            "mass": self.mass,
            "position": self.position,
            "velocity": self.velocity,
            "member_ids": self.member_ids,
        }
        if ptype == "halos" and field in data:
            return data[field]

        raise YTFieldNotFound((ptype, field), dataset=self.ds)

    @cached_property
    def ihalo(self):
        halo_id = self.particle_identifier
        halo_ids = list(self.halo_ds.r["halos", "particle_identifier"].astype("i8"))
        ihalo = halo_ids.index(halo_id)

        assert halo_ids[ihalo] == halo_id

        return ihalo

    @property
    def mass(self):
        return self.halo_ds.r["halos", "particle_mass"][self.ihalo]

    @property
    def position(self):
        return self.halo_ds.r["halos", "particle_position"][self.ihalo]

    @property
    def velocity(self):
        return self.halo_ds.r["halos", "particle_velocity"][self.ihalo]

    @property
    def member_ids(self):
        return self.halo_ds.index.get_member(self.particle_identifier)
