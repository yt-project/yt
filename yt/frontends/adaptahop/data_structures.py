"""
Data structures for AdaptaHOP frontend.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import stat
import os
import re

from .fields import \
    AdaptaHOPFieldInfo

from yt.data_objects.static_output import \
    Dataset
from yt.frontends.halo_catalog.data_structures import \
    HaloCatalogFile
from yt.funcs import \
    setdefaultattr
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.utilities.cython_fortran_utils import FortranFile
from yt.units import Mpc

from .definitions import \
    HEADER_ATTRIBUTES

class AdaptaHOPBinaryFile(HaloCatalogFile):
    def __init__(self, ds, io, filename, file_id):
        super(AdaptaHOPBinaryFile, self).__init__(ds, io, filename, file_id)

    def _read_particle_positions(self, ptype, f=None):
        raise NotImplementedError

class AdaptaHOPDataset(Dataset):
    _index_class = ParticleIndex
    _file_class = AdaptaHOPBinaryFile
    _field_info_class = AdaptaHOPFieldInfo

    _code_length_to_Mpc = Mpc.to('cm').value / 3.08e24

    def __init__(self, filename, dataset_type="adaptahop_binary",
                 n_ref = 16, over_refine_factor = 1,
                 units_override=None, unit_system="cgs",
                 parent_ds=None
                 ):
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        self.parent_ds = parent_ds
        super(AdaptaHOPDataset, self).__init__(filename, dataset_type,
                                              units_override=units_override,
                                              unit_system=unit_system)

    def _set_code_unit_attributes(self):
        setdefaultattr(self, 'length_unit', self.quan(self._code_length_to_Mpc, "Mpc"))
        setdefaultattr(self, 'mass_unit', self.quan(1e11, "Msun"))
        setdefaultattr(self, 'velocity_unit', self.quan(1.0, "km / s"))
        setdefaultattr(self, 'time_unit', self.length_unit / self.velocity_unit)

    def _parse_parameter_file(self):
        with FortranFile(self.parameter_filename) as fpu:
            params = fpu.read_attrs(HEADER_ATTRIBUTES)
        self.dimensionality = 3
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        # Domain related things
        self.filename_template = self.parameter_filename
        self.file_count = 1
        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(3, "int32") * nz

        # Set things up
        self.cosmological_simulation = 1
        self.current_redshift = (1.0 / params['aexp']) - 1.0
        self.omega_matter = params['omega_t']
        self.current_time = self.quan(params['age'], 'Gyr')
        self.omega_lambda = 0.724  # hard coded if not inferred from parent ds
        self.hubble_constant = 0.7  # hard coded if not inferred from parent ds
        self.periodicity = (True, True, True)
        self.particle_types = ("halos")
        self.particle_types_raw = ("halos")

        # Inherit stuff from parent ds
        if not self.parent_ds:
            raise Exception

        for k in ('omega_lambda', 'hubble_constant', 'omega_matter', 'omega_radiation'):
            setattr(self, k, getattr(self.parent_ds, k))

        self.domain_left_edge = np.array([0., 0., 0.])
        self.domain_right_edge = self.parent_ds.domain_right_edge.to('Mpc').value * \
            self._code_length_to_Mpc

        self.parameters.update(params)

    @classmethod
    def _is_valid(self, *args, **kwargs):
        dirname, fname = os.path.split(args[0])
        if not fname.startswith('tree_bricks') or not re.match('^tree_bricks\d{3}$', fname):
            return False
        return True

    def get_halo(self, halo_id, ptype='DM'):
        """
        Returns the smallest sphere that contains all the particle of the halo.
        """
        parent_ds = self.parent_ds
        ad = self.all_data()
        halo_ids = ad['halos', 'particle_identifier'].astype(int)
        ihalo = np.searchsorted(halo_ids, halo_id)
        assert halo_ids[ihalo] == halo_id

        halo_pos = ad['halos', 'particle_position'][ihalo, :].to('Mpc').value
        halo_vel = ad['halos', 'particle_velocity'][ihalo, :]
        halo_radius = ad['halos', 'r'][ihalo].to('Mpc').value

        members = self.index.io.members(ihalo)
        ok = False
        f = 1/1.1
        while not ok:
            f *= 1.1
            sph = parent_ds.sphere(
                parent_ds.arr(halo_pos, 'Mpc'),
                parent_ds.arr(f * halo_radius, 'Mpc'))

            part_ids = sph[ptype, 'particle_identity'].astype(int)

            ok = len(np.lib.arraysetops.setdiff1d(members, part_ids)) == 0

        # Set bulk velocity
        sph.set_field_parameter('bulk_velocity', (halo_vel.to('km/s').value, 'km/s'))

        # Build subregion that only contains halo particles
        reg = sph.cut_region(
            ['np.in1d(obj["io", "particle_identity"].astype(int), members)'],
            locals=dict(members=members, np=np))
        return (members, sph, reg)
