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
from numbers import \
    Number as numeric_type
import numpy as np
import os
import stat
import time
import weakref

from .fields import \
    YTDataContainerFieldInfo, \
    YTGridFieldInfo

from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.data_objects.static_output import \
    Dataset, \
    ParticleFile
from yt.extern.six import \
    iteritems, \
    string_types
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.utilities.cosmology import Cosmology
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

class YTDataContainerDataset(Dataset):
    _index_class = ParticleIndex
    _file_class = YTDataHDF5File
    _field_info_class = YTDataContainerFieldInfo
    _suffix = ".h5"

    def __init__(self, filename, dataset_type="ytdatacontainer_hdf5",
                 n_ref = 16, over_refine_factor = 1, units_override=None):
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        super(YTDataContainerDataset, self).__init__(filename, dataset_type,
                                            units_override=units_override)

    def _parse_parameter_file(self):
        with h5py.File(self.parameter_filename, "r") as f:
            hvals = dict((key, f.attrs[key]) for key in f.attrs.keys())
            self.particle_types_raw = tuple(f.keys())
        self.particle_types = self.particle_types_raw
        self.dimensionality = 3
        self.refine_by = 2
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        prefix = ".".join(self.parameter_filename.rsplit(".", 2)[:-2])
        self.filename_template = self.parameter_filename
        self.file_count = 1
        for attr in ["cosmological_simulation", "current_time", "current_redshift",
                     "hubble_constant", "omega_matter", "omega_lambda",
                     "domain_left_edge", "domain_right_edge"]:
            setattr(self, attr, hvals[attr])
        self.periodicity = (True, True, True)

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

class YTGrid(AMRGridPatch):
    pass

class YTGridHierarchy(GridIndex):
    grid = YTGrid

    def _count_grids(self):
        self.num_grids = 1

class YTGridDataset(Dataset):
    _index_class = YTGridHierarchy
    _field_info_class = YTGridFieldInfo
    _dataset_type = 'ytgridhdf5'
    geometry = "cartesian"

    def __init__(self, filename):
        Dataset.__init__(self, filename, self._dataset_type)

    def _parse_parameter_file(self):
        self.dimensionality = 3
        self.refine_by = 2
        self.unique_identifier = time.time()
        with h5py.File(self.parameter_filename, "r") as f:
            for attr, value in f.attrs.items():
                setattr(self, attr, value)

        # correct domain dimensions for the covering grid dimension
        self.base_domain_left_edge = self.domain_left_edge
        self.base_domain_right_edge = self.domain_right_edge
        self.base_domain_dimensions = self.domain_dimensions
        dx = (self.domain_right_edge - self.domain_left_edge) / \
          (self.domain_dimensions * self.refine_by**self.level)
        self.domain_left_edge = self.left_edge
        self.domain_right_edge = self.domain_left_edge + \
          self.ActiveDimensions * dx
        self.domain_dimensions = self.ActiveDimensions

    def __repr__(self):
        return "ytGrid: %s" % self.parameter_filename

    def _set_code_unit_attributes(self):
        attrs = ('length_unit', 'mass_unit', 'time_unit',
                 'velocity_unit', 'magnetic_unit')
        cgs_units = ('cm', 'g', 's', 'cm/s', 'gauss')
        base_units = np.ones(len(attrs))
        for unit, attr, cgs_unit in zip(base_units, attrs, cgs_units):
            if isinstance(unit, string_types):
                uq = self.quan(1.0, unit)
            elif isinstance(unit, numeric_type):
                uq = self.quan(unit, cgs_unit)
            elif isinstance(unit, YTQuantity):
                uq = unit
            elif isinstance(unit, tuple):
                uq = self.quan(unit[0], unit[1])
            else:
                raise RuntimeError("%s (%s) is invalid." % (attr, unit))
            setattr(self, attr, uq)

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not args[0].endswith(".h5"): return False
        with h5py.File(args[0], "r") as f:
            if "data_type" in f.attrs and \
              f.attrs["data_type"] in ["yt_grid_data"]:
                return True
        return False
