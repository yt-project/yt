"""
Data structures for Rockstar frontend.




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
import glob
import os

from .fields import \
    RockstarFieldInfo

from yt.data_objects.static_output import \
    Dataset, \
    ParticleFile
from yt.funcs import \
    setdefaultattr
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.utilities.cosmology import Cosmology
import yt.utilities.fortran_utils as fpu

from .definitions import \
    header_dt

class RockstarBinaryFile(ParticleFile):
    def __init__(self, ds, io, filename, file_id):
        with open(filename, "rb") as f:
            self.header = fpu.read_cattrs(f, header_dt, "=")
            self._position_offset = f.tell()
            f.seek(0, os.SEEK_END)
            self._file_size = f.tell()

        super(RockstarBinaryFile, self).__init__(ds, io, filename, file_id)

class RockstarDataset(Dataset):
    _index_class = ParticleIndex
    _file_class = RockstarBinaryFile
    _field_info_class = RockstarFieldInfo
    _suffix = ".bin"

    def __init__(self, filename, dataset_type="rockstar_binary",
                 n_ref = 16, over_refine_factor = 1,
                 units_override=None, unit_system="cgs"):
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        super(RockstarDataset, self).__init__(filename, dataset_type,
                                              units_override=units_override,
                                              unit_system=unit_system)

    def _parse_parameter_file(self):
        with open(self.parameter_filename, "rb") as f:
            hvals = fpu.read_cattrs(f, header_dt)
            hvals.pop("unused")
        self.dimensionality = 3
        self.refine_by = 2
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        prefix = ".".join(self.parameter_filename.rsplit(".", 2)[:-2])
        self.filename_template = "%s.%%(num)s%s" % (prefix, self._suffix)
        self.file_count = len(glob.glob(prefix + ".*" + self._suffix))
        
        # Now we can set up things we already know.
        self.cosmological_simulation = 1
        self.current_redshift = (1.0 / hvals['scale']) - 1.0
        self.hubble_constant = hvals['h0']
        self.omega_lambda = hvals['Ol']
        self.omega_matter = hvals['Om']
        cosmo = Cosmology(self.hubble_constant,
                          self.omega_matter, self.omega_lambda)
        self.current_time = cosmo.hubble_time(self.current_redshift).in_units("s")
        self.periodicity = (True, True, True)
        self.particle_types = ("halos")
        self.particle_types_raw = ("halos")

        self.domain_left_edge = np.array([0.0,0.0,0.0])
        self.domain_right_edge = np.array([hvals['box_size']] * 3)

        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(3, "int32") * nz
        self.parameters.update(hvals)

    def _set_code_unit_attributes(self):
        z = self.current_redshift
        setdefaultattr(self, 'length_unit', self.quan(1.0 / (1.0+z), "Mpc / h"))
        setdefaultattr(self, 'mass_unit', self.quan(1.0, "Msun / h"))
        setdefaultattr(self, 'velocity_unit', self.quan(1.0, "km / s"))
        setdefaultattr(self, 'time_unit', self.length_unit / self.velocity_unit)

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not args[0].endswith(".bin"): return False
        with open(args[0], "rb") as f:
            header = fpu.read_cattrs(f, header_dt)
            if header['magic'] == 18077126535843729616:
                return True
        return False
