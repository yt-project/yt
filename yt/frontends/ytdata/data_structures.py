"""
Data structures for YTData frontend.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
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

from .fields import \
    YTDataFieldInfo

from yt.utilities.cosmology import Cosmology
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.data_objects.static_output import \
    Dataset, \
    ParticleFile
import yt.utilities.fortran_utils as fpu
from yt.units.yt_array import \
    YTArray, \
    YTQuantity
    
class YTDataHDF5File(ParticleFile):
    def __init__(self, ds, io, filename, file_id):
        with h5py.File(filename, "r") as f:
            self.header = dict((field, f.attrs[field]) \
                               for field in f.attrs.keys())

        super(YTDataHDF5File, self).__init__(ds, io, filename, file_id)
    
class YTDataDataset(Dataset):
    _index_class = ParticleIndex
    _file_class = YTDataHDF5File
    _field_info_class = YTDataFieldInfo
    _suffix = ".h5"

    def __init__(self, filename, dataset_type="ytdata_hdf5",
                 n_ref = 16, over_refine_factor = 1, units_override=None):
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        super(YTDataDataset, self).__init__(filename, dataset_type,
                                                 units_override=units_override)

    def _parse_parameter_file(self):
        with h5py.File(self.parameter_filename, "r") as f:
            hvals = dict((key, f.attrs[key]) for key in f.attrs.keys())
        self.dimensionality = 3
        self.refine_by = 2
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        prefix = ".".join(self.parameter_filename.rsplit(".", 2)[:-2])
        self.filename_template = "%s.%%(num)s%s" % (prefix, self._suffix)
        self.file_count = len(glob.glob(prefix + "*" + self._suffix))

        for attr in ["cosmological_simulation", "current_time", "current_redshift",
                     "hubble_constant", "omega_matter", "omega_lambda",
                     "domain_left_edge", "domain_right_edge"]:
            setattr(self, attr, hvals[attr])
        self.periodicity = (True, True, True)
        self.particle_types = ("grid")
        self.particle_types_raw = ("grid")

        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(3, "int32") * nz
        self.parameters.update(hvals)

    def _set_code_unit_attributes(self):
        self.length_unit = self.quan(1.0, "cm")
        self.mass_unit = self.quan(1.0, "g")
        self.velocity_unit = self.quan(1.0, "cm / s")
        self.time_unit = self.quan(1.0, "s")

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not args[0].endswith(".h5"): return False
        with h5py.File(args[0], "r") as f:
            if "data_type" in f.attrs and \
              f.attrs["data_type"] in ["light_ray",
                                       "yt_array_data",
                                       "yt_data_container"]:
                return True
        return False
