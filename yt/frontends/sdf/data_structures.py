import os
import sys

import numpy as np

from yt.data_objects.static_output import ParticleDataset, ParticleFile
from yt.funcs import setdefaultattr
from yt.geometry.particle_geometry_handler import ParticleIndex
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _requests as requests
from yt.utilities.sdf import HTTPSDFRead, SDFIndex, SDFRead

from .fields import SDFFieldInfo

if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    from yt._maintenance.backports import cached_property
# currently specified by units_2HOT == 2 in header
# in future will read directly from file
units_2HOT_v2_length = 3.08567802e21
units_2HOT_v2_mass = 1.98892e43
units_2HOT_v2_time = 3.1558149984e16


class SDFFile(ParticleFile):
    pass


class SDFDataset(ParticleDataset):
    _index_class = ParticleIndex
    _file_class = SDFFile
    _field_info_class = SDFFieldInfo
    _particle_mass_name = None
    _particle_coordinates_name = None
    _particle_velocity_name = None
    _skip_cache = True
    _subspace = False

    def __init__(
        self,
        filename,
        dataset_type="sdf_particles",
        index_order=None,
        index_filename=None,
        bounding_box=None,
        sdf_header=None,
        midx_filename=None,
        midx_header=None,
        midx_level=None,
        field_map=None,
        units_override=None,
        unit_system="cgs",
    ):
        if bounding_box is not None:
            # This ensures that we know a bounding box has been applied
            self._domain_override = True
            self._subspace = True
            bbox = np.array(bounding_box, dtype="float64")
            if bbox.shape == (2, 3):
                bbox = bbox.transpose()
            self.domain_left_edge = bbox[:, 0]
            self.domain_right_edge = bbox[:, 1]
        else:
            self.domain_left_edge = self.domain_right_edge = None
        self.sdf_header = sdf_header
        self.midx_filename = midx_filename
        self.midx_header = midx_header
        self.midx_level = midx_level
        if field_map is None:
            field_map = {}
        self._field_map = field_map
        prefix = ""
        if self.midx_filename is not None:
            prefix += "midx_"
        if filename.startswith("http"):
            prefix += "http_"
        dataset_type = prefix + "sdf_particles"
        super().__init__(
            filename,
            dataset_type=dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            index_order=index_order,
            index_filename=index_filename,
        )

    def _parse_parameter_file(self):
        if self.parameter_filename.startswith("http"):
            sdf_class = HTTPSDFRead
        else:
            sdf_class = SDFRead
        self.sdf_container = sdf_class(self.parameter_filename, header=self.sdf_header)

        # Reference
        self.parameters = self.sdf_container.parameters
        self.dimensionality = 3
        self.refine_by = 2

        if self.domain_left_edge is None or self.domain_right_edge is None:
            R0 = self.parameters["R0"]
            if "offset_center" in self.parameters and self.parameters["offset_center"]:
                self.domain_left_edge = np.array([0, 0, 0], dtype=np.float64)
                self.domain_right_edge = np.array(
                    [2.0 * self.parameters.get(f"R{ax}", R0) for ax in "xyz"],
                    dtype=np.float64,
                )
            else:
                self.domain_left_edge = np.array(
                    [-self.parameters.get(f"R{ax}", R0) for ax in "xyz"],
                    dtype=np.float64,
                )
                self.domain_right_edge = np.array(
                    [+self.parameters.get(f"R{ax}", R0) for ax in "xyz"],
                    dtype=np.float64,
                )
            self.domain_left_edge *= self.parameters.get("a", 1.0)
            self.domain_right_edge *= self.parameters.get("a", 1.0)

        self.domain_dimensions = np.ones(3, "int32")
        if self.parameters.get("do_periodic", False):
            self._periodicity = (True, True, True)
        else:
            self._periodicity = (False, False, False)

        self.cosmological_simulation = 1

        self.current_redshift = self.parameters.get("redshift", 0.0)
        self.omega_lambda = self.parameters["Omega0_lambda"]
        self.omega_matter = self.parameters["Omega0_m"]
        if "Omega0_fld" in self.parameters:
            self.omega_lambda += self.parameters["Omega0_fld"]
        if "Omega0_r" in self.parameters:
            # not correct, but most codes can't handle Omega0_r
            self.omega_matter += self.parameters["Omega0_r"]
        self.hubble_constant = self.parameters["h_100"]
        self.current_time = units_2HOT_v2_time * self.parameters.get("tpos", 0.0)
        mylog.info("Calculating time to be %0.3e seconds", self.current_time)
        self.filename_template = self.parameter_filename
        self.file_count = 1

    @cached_property
    def midx(self):
        if self.midx_filename is None:
            raise RuntimeError("SDF index0 file not supplied in load.")

        if "http" in self.midx_filename:
            sdf_class = HTTPSDFRead
        else:
            sdf_class = SDFRead
        indexdata = sdf_class(self.midx_filename, header=self.midx_header)
        return SDFIndex(self.sdf_container, indexdata, level=self.midx_level)

    def _set_code_unit_attributes(self):
        setdefaultattr(
            self,
            "length_unit",
            self.quan(1.0, self.parameters.get("length_unit", "kpc")),
        )
        setdefaultattr(
            self,
            "velocity_unit",
            self.quan(1.0, self.parameters.get("velocity_unit", "kpc/Gyr")),
        )
        setdefaultattr(
            self, "time_unit", self.quan(1.0, self.parameters.get("time_unit", "Gyr"))
        )
        mass_unit = self.parameters.get("mass_unit", "1e10 Msun")
        if " " in mass_unit:
            factor, unit = mass_unit.split(" ")
        else:
            factor = 1.0
            unit = mass_unit
        setdefaultattr(self, "mass_unit", self.quan(float(factor), unit))

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        sdf_header = kwargs.get("sdf_header", filename)
        if sdf_header.startswith("http"):
            try:
                hreq = requests.get(sdf_header, stream=True)
            except ImportError:
                # requests is not installed
                return False
            if hreq.status_code != 200:
                return False
            # Grab a whole 4k page.
            line = next(hreq.iter_content(4096))
        elif os.path.isfile(sdf_header):
            try:
                with open(sdf_header, encoding="ISO-8859-1") as f:
                    line = f.read(10).strip()
            except PermissionError:
                return False
        else:
            return False
        return line.startswith("# SDF")
