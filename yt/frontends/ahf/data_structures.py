"""
AHF data structures



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import glob
import os
import stat

import numpy as np

from yt.data_objects.static_output import \
    Dataset, \
    ParticleFile
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.utilities.cosmology import \
    Cosmology

from .fields import AHFHalosFieldInfo


class AHFHalosFile(ParticleFile):
    def __init__(self, ds, io, filename, file_id):
        root, _ = os.path.splitext(filename)
        candidates = glob.glob(root + '*.AHF_halos')
        if len(candidates) == 1:
            filename = candidates[0]
        else:
            raise ValueError('Too many AHF_halos files.')
        names = self._read_column_names(filename)
        self.data = np.genfromtxt(filename, names=names)
        super(AHFHalosFile, self).__init__(ds, io, filename, file_id)

    def _read_column_names(self, filename):
        with open(filename) as f:
            line = f.readline()
            # Remove leading '#'
            line = line[1:]
            names = line.split()
            # Remove trailing '()'
            names = [name.split('(')[0] for name in names]
            return names


class AHFHalosDataset(Dataset):
    _index_class = ParticleIndex
    _file_class = AHFHalosFile
    _field_info_class = AHFHalosFieldInfo

    def __init__(self, filename, dataset_type='ahf',
                 n_ref=16, over_refine_factor=1, units_override=None,
                 unit_system='cgs'):
        root, _ = os.path.splitext(filename)
        self.log_filename = root + '.log'

        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        super(AHFHalosDataset, self).__init__(
            filename, dataset_type=dataset_type,
            units_override=units_override, unit_system=unit_system
        )

    def _set_code_unit_attributes(self):
        self.length_unit = self.quan(1.0, 'kpccm/h')
        self.mass_unit = self.quan(1.0, 'Msun/h')
        self.time_unit = self.quan(1.0, '')
        self.velocity_unit = self.quan(1.0, 'km/s')

    def _parse_parameter_file(self):
        # This needs to set up the following items.  Note that these are all
        # assumed to be in code units; domain_left_edge and domain_right_edge
        # will be converted to YTArray automatically at a later time.
        # This includes the cosmological parameters.
        #
        #   self.unique_identifier      <= unique identifier for the dataset
        #                                  being read (e.g., UUID or ST_CTIME)
        #   self.parameters             <= full of code-specific items of use
        #   self.domain_left_edge       <= array of float64
        #   self.domain_right_edge      <= array of float64
        #   self.dimensionality         <= int
        #   self.domain_dimensions      <= array of int64
        #   self.periodicity            <= three-element tuple of booleans
        #   self.current_time           <= simulation time in code units
        #
        # We also set up cosmological information.  Set these to zero if
        # non-cosmological.
        #
        #   self.cosmological_simulation    <= int, 0 or 1
        #   self.current_redshift           <= float
        #   self.omega_lambda               <= float
        #   self.omega_matter               <= float
        #   self.hubble_constant            <= float

        # Read all parameters.
        simu = self._read_log_simu()
        param = self._read_parameter()

        # Set up general information.
        self.filename_template = self.parameter_filename
        self.file_count = 1
        self.parameters.update(param)
        self.particle_types = ('halos')
        self.particle_types_raw = ('halos')
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # Set up geometrical information.
        self.refine_by = 2
        self.dimensionality = 3
        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(self.dimensionality, "int32") * nz
        self.domain_left_edge = np.array([0.0, 0.0, 0.0])
        self.domain_right_edge = np.array([simu['boxsize']] * 3)
        self.periodicity = (True, True, True)

        # Set up cosmological information.
        self.cosmological_simulation = 1
        self.current_redshift = param['z']
        self.hubble_constant = param['Hubble(z)']
        self.omega_lambda = simu['lambda0']
        self.omega_matter = simu['omega0']
        cosmo = Cosmology(self.hubble_constant,
                          self.omega_matter, self.omega_lambda)
        self.current_time = cosmo.hubble_time(param['z']).in_units('s')

    @classmethod
    def _is_valid(self, *args, **kwargs):
        filename = args[0]
        if not filename.endswith('.parameter'):
            return False
        with open(filename, 'r') as f:
            if f.readlines()[11].startswith('AHF'):
                return True
        return False

    # Helper methods

    def _read_log_simu(self):
        simu = {}
        with open(self.log_filename) as f:
            for l in f:
                if l.startswith('simu.'):
                    name, val = l.split(':')
                    key = name.strip().split('.')[1]
                    try:
                        val = float(val)
                    except:
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
                    except:
                        pass
        return param
