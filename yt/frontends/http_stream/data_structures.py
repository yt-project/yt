import json
import sys
import time

import numpy as np

from yt.data_objects.static_output import ParticleDataset, ParticleFile
from yt.frontends.sph.fields import SPHFieldInfo
from yt.geometry.particle_geometry_handler import ParticleIndex
from yt.utilities.on_demand_imports import _requests as requests

if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    from yt._maintenance.backports import cached_property


class HTTPParticleFile(ParticleFile):
    pass


class HTTPStreamDataset(ParticleDataset):
    _index_class = ParticleIndex
    _file_class = HTTPParticleFile
    _field_info_class = SPHFieldInfo
    _particle_mass_name = "Mass"
    _particle_coordinates_name = "Coordinates"
    _particle_velocity_name = "Velocities"
    filename_template = ""

    def __init__(
        self,
        base_url,
        dataset_type="http_particle_stream",
        unit_system="cgs",
        index_order=None,
        index_filename=None,
    ):
        self.base_url = base_url
        super().__init__(
            "",
            dataset_type=dataset_type,
            unit_system=unit_system,
            index_order=index_order,
            index_filename=index_filename,
        )

    def __str__(self):
        return self.base_url

    @cached_property
    def unique_identifier(self) -> str:
        return str(self.parameters.get("unique_identifier", time.time()))

    def _parse_parameter_file(self):
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"

        # Here's where we're going to grab the JSON index file
        hreq = requests.get(self.base_url + "/yt_index.json")
        if hreq.status_code != 200:
            raise RuntimeError
        header = json.loads(hreq.content)
        header["particle_count"] = {
            int(k): header["particle_count"][k] for k in header["particle_count"]
        }
        self.parameters = header

        # Now we get what we need
        self.domain_left_edge = np.array(header["domain_left_edge"], "float64")
        self.domain_right_edge = np.array(header["domain_right_edge"], "float64")
        self.domain_dimensions = np.ones(3, "int32")
        self._periodicity = (True, True, True)

        self.current_time = header["current_time"]
        self.cosmological_simulation = int(header["cosmological_simulation"])
        for attr in (
            "current_redshift",
            "omega_lambda",
            "omega_matter",
            "hubble_constant",
        ):
            setattr(self, attr, float(header[attr]))

        self.file_count = header["num_files"]

    def _set_units(self):
        length_unit = float(self.parameters["units"]["length"])
        time_unit = float(self.parameters["units"]["time"])
        mass_unit = float(self.parameters["units"]["mass"])
        density_unit = mass_unit / length_unit**3
        velocity_unit = length_unit / time_unit
        self._unit_base = {}
        self._unit_base["cm"] = 1.0 / length_unit
        self._unit_base["s"] = 1.0 / time_unit
        super()._set_units()
        self.conversion_factors["velocity"] = velocity_unit
        self.conversion_factors["mass"] = mass_unit
        self.conversion_factors["density"] = density_unit

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if not filename.startswith("http://"):
            return False
        try:
            return requests.get(filename + "/yt_index.json").status_code == 200
        except ImportError:
            # requests is not installed
            return False
