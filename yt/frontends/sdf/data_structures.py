"""
Data structures for a generic SDF frontend




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import numpy as np
import stat
import weakref
import struct
import glob
import time
import os
import types

from yt.utilities.logger import ytLogger as mylog
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.data_objects.static_output import \
    Dataset, ParticleFile
from yt.utilities.physical_constants import \
    G, \
    cm_per_kpc, \
    mass_sun_cgs
from yt.utilities.cosmology import Cosmology
from .fields import \
    SDFFieldInfo
from .io import \
    IOHandlerSDF, \
    SDFRead

class SDFFile(ParticleFile):
    def __init__(self, pf, io, filename, file_id):
        pass

class SDFDataset(Dataset):
    _index_class = ParticleIndex
    _file_class = SDFFile
    _field_info_class = SDFFieldInfo
    _particle_mass_name = None
    _particle_coordinates_name = None
    _particle_velocity_name = None

    def __init__(self, filename, dataset_type = "sdf_particles",
                 n_ref = 64, over_refine_factor = 1,
                 bounding_box = None):
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        if bounding_box is not None:
            bbox = np.array(bounding_box, dtype="float64")
            if bbox.shape == (2, 3):
                bbox = bbox.transpose()
            self.domain_left_edge = bbox[:,0]
            self.domain_right_edge = bbox[:,1]
        else:
            self.domain_left_edge = self.domain_right_edge = None
        super(SDFDataset, self).__init__(filename, dataset_type)

    def _parse_parameter_file(self):
        self.sdf_container = SDFRead(self.parameter_filename)
        # Reference
        self.parameters = self.sdf_container.parameters
        self.dimensionality = 3
        self.refine_by = 2
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        if None in (self.domain_left_edge, self.domain_right_edge):
            self.domain_left_edge = np.array([
              -self.parameters["R%s" % ax] for ax in 'xyz'])
            self.domain_right_edge = np.array([
              +self.parameters["R%s" % ax] for ax in 'xyz'])
            self.domain_left_edge *= self.parameters["a"]
            self.domain_right_edge *= self.parameters["a"]

        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(3, "int32") * nz
        self.periodicity = (True, True, True)

        self.cosmological_simulation = 1

        self.current_redshift = self.parameters["redshift"]
        self.omega_lambda = self.parameters["Omega0_lambda"]
        self.omega_matter = self.parameters["Omega0_m"]
        self.hubble_constant = self.parameters["hubble"]
        # Now we calculate our time based on the cosmology.
        cosmo = Cosmology(self.hubble_constant,
                          self.omega_matter, self.omega_lambda)
        self.current_time = cosmo.hubble_time(self.current_redshift)
        mylog.info("Calculating time to be %0.3e seconds", self.current_time)

    def _set_code_unit_attributes(self):
        self.length_unit = self.quan(1.0, "kpc")
        self.velocity_unit = self.quan(1.0, "kpc/Gyr")
        self.time_unit = self.quan(1.0, "Gyr")
        self.mass_unit = self.quan(1e10, "Msun")

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        with open(args[0], "r") as f:
            line = f.readline().strip()
            return line == "# SDF 1.0"
