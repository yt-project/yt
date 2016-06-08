"""
Data structures for HTTPStream frontend.




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import json
import numpy as np
import time

from yt.data_objects.static_output import \
    ParticleFile
from yt.frontends.sph.data_structures import \
    ParticleDataset
from yt.frontends.sph.fields import \
    SPHFieldInfo
from yt.funcs import \
    get_requests
from yt.geometry.particle_geometry_handler import \
    ParticleIndex

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
    
    def __init__(self, base_url,
                 dataset_type = "http_particle_stream",
                 n_ref = 64, over_refine_factor=1, 
                 unit_system="cgs"):
        if get_requests() is None:
            raise ImportError(
                "This functionality depends on the requests package")
        self.base_url = base_url
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        super(HTTPStreamDataset, self).__init__("", dataset_type, 
                                                unit_system=unit_system)

    def __repr__(self):
        return self.base_url

    def _parse_parameter_file(self):
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"

        # Here's where we're going to grab the JSON index file
        requests = get_requests()
        hreq = requests.get(self.base_url + "/yt_index.json")
        if hreq.status_code != 200:
            raise RuntimeError
        header = json.loads(hreq.content)
        header['particle_count'] = dict((int(k), header['particle_count'][k])
            for k in header['particle_count'])
        self.parameters = header

        # Now we get what we need
        self.domain_left_edge = np.array(header['domain_left_edge'], "float64")
        self.domain_right_edge = np.array(header['domain_right_edge'], "float64")
        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(3, "int32") * nz
        self.periodicity = (True, True, True)

        self.current_time = header['current_time']
        self.unique_identifier = header.get("unique_identifier", time.time())
        self.cosmological_simulation = int(header['cosmological_simulation'])
        for attr in ('current_redshift', 'omega_lambda', 'omega_matter',
                     'hubble_constant'):
            setattr(self, attr, float(header[attr]))

        self.file_count = header['num_files']

    def _set_units(self):
        length_unit = float(self.parameters['units']['length'])
        time_unit = float(self.parameters['units']['time'])
        mass_unit = float(self.parameters['units']['mass'])
        density_unit = mass_unit / length_unit ** 3
        velocity_unit = length_unit / time_unit
        self._unit_base = {}
        self._unit_base['cm'] = 1.0/length_unit
        self._unit_base['s'] = 1.0/time_unit
        super(HTTPStreamDataset, self)._set_units()
        self.conversion_factors["velocity"] = velocity_unit
        self.conversion_factors["mass"] = mass_unit
        self.conversion_factors["density"] = density_unit

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not args[0].startswith("http://"):
            return False
        requests = get_requests()
        if requests is None:
            return False
        hreq = requests.get(args[0] + "/yt_index.json")
        if hreq.status_code == 200:
            return True
        return False
