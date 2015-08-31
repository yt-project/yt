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
    YTProjectionFieldInfo, \
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
from yt.utilities.logger import \
    ytLogger as mylog
from yt.utilities.cosmology import Cosmology
import yt.utilities.fortran_utils as fpu
from yt.units.yt_array import \
    YTArray, \
    YTQuantity

_grid_data_containers = ["abritrary_grid",
                         "covering_grid",
                         "smoothed_covering_grid"]

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
        self.refine_by = 2
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        prefix = ".".join(self.parameter_filename.rsplit(".", 2)[:-2])
        self.filename_template = self.parameter_filename
        self.file_count = 1
        for attr in ["cosmological_simulation", "current_time", "current_redshift",
                     "hubble_constant", "omega_matter", "omega_lambda",
                     "dimensionality", "domain_left_edge", "domain_right_edge"]:
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
            data_type = f.attrs.get("data_type", None)
            if data_type is None:
                return False
            if data_type in ["yt_light_ray", "yt_array_data"]:
                return True
            if data_type == "yt_data_container" and \
              f.attrs.get("container_type", None) not in \
              _grid_data_containers:
                return True
        return False

class YTProjectionDataset(YTDataContainerDataset):
    _field_info_class = YTProjectionFieldInfo

    def __init__(self, *args, **kwargs):
        super(YTProjectionDataset, self).__init__(
            *args, dataset_type="ytprojection_hdf5", **kwargs)

    def _parse_parameter_file(self):
        super(YTProjectionDataset, self)._parse_parameter_file()
        self.axis = self.parameters["axis"]
        self.weight_field = self.parameters["weight_field"]
        if isinstance(self.weight_field, str) and \
          self.weight_field == "None":
            self.weight_field = None
        elif isinstance(self.weight_field, np.ndarray):
            self.weight_field = tuple(self.weight_field)

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not args[0].endswith(".h5"): return False
        with h5py.File(args[0], "r") as f:
            data_type = f.attrs.get("data_type", None)
            if data_type == "yt_data_container" and \
              f.attrs.get("container_type", None) == "proj":
                return True
        return False

class YTGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, id, index):
        AMRGridPatch.__init__(self, id, filename=None, index=index)
        self._children_ids = []
        self._parent_id = -1
        self.Level = 0
        self.LeftEdge = self.index.ds.domain_left_edge
        self.RightEdge = self.index.ds.domain_right_edge

    @property
    def Parent(self):
        return None

    @property
    def Children(self):
        return []

class YTGridHierarchy(GridIndex):
    grid = YTGrid

    def __init__(self, ds, dataset_type = None):
        self.dataset_type = dataset_type
        self.float_type = 'float64'
        self.dataset = weakref.proxy(ds)
        self.directory = os.getcwd()
        GridIndex.__init__(self, ds, dataset_type)

    def _count_grids(self):
        self.num_grids = 1

    def _parse_index(self):
        self.grid_dimensions[:] = self.ds.domain_dimensions
        self.grid_left_edge[:] = self.ds.domain_left_edge
        self.grid_right_edge[:] = self.ds.domain_right_edge
        self.grid_levels[:] = np.zeros(self.num_grids)
        self.grid_procs = np.zeros(self.num_grids)
        self.grid_particle_count[:] = sum(self.ds.num_particles.values())
        self.grids = []
        for id in range(self.num_grids):
            self.grids.append(self.grid(id, self))
            self.grids[id].Level = self.grid_levels[id, 0]
        self.max_level = self.grid_levels.max()
        temp_grids = np.empty(self.num_grids, dtype='object')
        for i, grid in enumerate(self.grids):
            grid.filename = self.ds.parameter_filename
            grid._prepare_grid()
            grid.proc_num = self.grid_procs[i]
            temp_grids[i] = grid
        self.grids = temp_grids

    def _populate_grid_objects(self):
        for g in self.grids:
            g._setup_dx()
        self.max_level = self.grid_levels.max()

    def _detect_output_fields(self):
        self.field_list = []
        self.ds.field_units = self.ds.field_units or {}
        with h5py.File(self.ds.parameter_filename, "r") as f:
            for group in f:
                for field in f[group]:
                    field_name = (str(group), str(field))
                    self.field_list.append(field_name)
                    self.ds.field_units[field_name] = \
                      f[group][field].attrs["units"]

class YTGridDataset(Dataset):
    _index_class = YTGridHierarchy
    _field_info_class = YTGridFieldInfo
    _dataset_type = 'ytgridhdf5'
    geometry = "cartesian"
    default_fluid_type = "grid"
    fluid_types = ("grid", "gas", "deposit", "index")

    def __init__(self, filename):
        Dataset.__init__(self, filename, self._dataset_type)

    def _parse_parameter_file(self):
        self.refine_by = 2
        self.unique_identifier = time.time()
        with h5py.File(self.parameter_filename, "r") as f:
            for attr, value in f.attrs.items():
                setattr(self, attr, value)
            self.num_particles = \
              dict([(group, f[group].attrs["num_elements"])
                    for group in f if group != self.default_fluid_type])
        self.particle_types_raw = tuple(self.num_particles.keys())
        self.particle_types = self.particle_types_raw

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
        self.periodicity = \
          np.abs(self.domain_left_edge -
                 self.base_domain_left_edge) < 0.5 * dx
        self.periodicity &= \
        np.abs(self.domain_right_edge -
               self.base_domain_right_edge) < 0.5 * dx

    def __repr__(self):
        return "ytGrid: %s" % self.parameter_filename

    def create_field_info(self):
        self.field_info = self._field_info_class(self, self.field_list)
        for ftype, field in self.field_list:
            if ftype == self.default_fluid_type:
                self.field_info.alias(
                    ("gas", field),
                    (self.default_fluid_type, field))
        super(YTGridDataset, self).create_field_info()

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
            data_type = f.attrs.get("data_type", None)
            if data_type == "yt_data_container" and \
              f.attrs.get("container_type", None) in \
              _grid_data_containers:
                return True
        return False
