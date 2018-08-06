"""
CF Radial data structures



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import os
import stat
import weakref

from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from .fields import CFRadialFieldInfo

import xarray

class CFRadialGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level, dimensions):
        super(CFRadialGrid, self).__init__(
            id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level
        self.ActiveDimensions = dimensions

    def __repr__(self):
        return "CFRadialGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class CFRadialHierarchy(GridIndex):
    grid = CFRadialGrid

    def __init__(self, ds, dataset_type='cf_radial'):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super(CFRadialHierarchy, self).__init__(
            ds, dataset_type)

    def _initialize_state_variables(self):
        super(CFRadialHierarchy, self)._initialize_state_variables()
        self.num_grids = 1

    def _detect_output_fields(self):
        # This needs to set a self.field_list that contains all the available,
        # on-disk fields. No derived fields should be defined here.
        # NOTE: Each should be a tuple, where the first element is the on-disk
        # fluid type or particle type.  Convention suggests that the on-disk
        # fluid type is usually the dataset_type and the on-disk particle type
        # (for a single population of particles) is "io".
        self.field_list = [
            ('cf_radial', 'reflectivity')
        ]

    def _count_grids(self):
        # This needs to set self.num_grids
        pass

    def _parse_index(self):
        # This needs to fill the following arrays, where N is self.num_grids:
        #   self.grid_left_edge         (N, 3) <= float64
        #   self.grid_right_edge        (N, 3) <= float64
        #   self.grid_dimensions        (N, 3) <= int
        #   self.grid_particle_count    (N, 1) <= int
        #   self.grid_levels            (N, 1) <= int
        #   self.grids                  (N, 1) <= grid objects
        #   self.max_level = self.grid_levels.max()
        self.grid_left_edge[0][:] = self.ds.domain_left_edge[:]
        self.grid_right_edge[0][:] = self.ds.domain_right_edge[:]
        self.grid_dimensions[0][:] = self.ds.domain_dimensions[:]
        self.grids = np.empty(self.num_grids, dtype='object')
        for i in range(self.num_grids):
            g = self.grid(i, self, self.grid_levels.flat[i],
                          self.grid_dimensions[i])
            g._prepare_grid()
            g._setup_dx()
            self.grids[i] = g
        self.max_level = 1
        pass

    def _populate_grid_objects(self):
        # For each grid, this must call:
        #   grid._prepare_grid()
        #   grid._setup_dx()
        # This must also set:
        #   grid.Children <= list of child grids
        #   grid.Parent   <= parent grid
        # This is handled by the frontend because often the children must be
        # identified.
        pass


class CFRadialDataset(Dataset):
    _index_class = CFRadialHierarchy
    _field_info_class = CFRadialFieldInfo

    def __init__(self, filename, dataset_type='cf_radial',
                 storage_filename=None,
                 units_override=None):
        self.fluid_types += ('cf_radial',)
        self._handle = xarray.open_dataset(filename)
        super(CFRadialDataset, self).__init__(
            filename, dataset_type,
            units_override=units_override)
        self.storage_filename = storage_filename
        # refinement factor between a grid and its subgrid
        # self.refine_by = 2

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        # self.length_unit = self.quan(1.0, "cm")
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "cm/s")
        # self.magnetic_unit = self.quan(1.0, "gauss")
        length_unit = self._handle.variables['x'].attrs['units']
        self.length_unit = self.quan(1.0, length_unit)
        self.mass_unit = self.quan(1.0, "kg")
        self.time_unit = self.quan(1.0, "s")


    def _parse_parameter_file(self):
        # This needs to set up the following items.  Note that these are all
        # assumed to be in code units; domain_left_edge and domain_right_edge
        # will be converted to YTArray automatically at a later time.
        # This includes the cosmological parameters.

        #
        #   self.unique_identifier      <= unique identifier for the dataset
        #                                  being read (e.g., UUID or ST_CTIME)
        self.unique_identifier = int(os.stat(
            self.parameter_filename)[stat.ST_CTIME])
        #   self.parameters             <= full of code-specific items of use
        self.parameters = {}

        coords = self._handle.coords
        x, y, z = [coords[d] for d in 'xyz']

        self.origin_latitude = self._handle.origin_latitude[0]
        self.origin_longitude = self._handle.origin_longitude[0]

        #   self.domain_left_edge       <= array of float64
        dle = [x.min(), y.min(), z.min()]
        self.domain_left_edge = np.array(dle)
        #   self.domain_right_edge      <= array of float64
        dre = [x.max(), y.max(), z.max()]
        self.domain_right_edge = np.array(dre)
        #   self.dimensionality         <= int
        self.dimensionality = 3
        #   self.domain_dimensions      <= array of int64
        dims = [self._handle.dims[d] for d in 'xyz']
        self.domain_dimensions = np.array(dims, dtype='int64')
        #   self.periodicity            <= three-element tuple of booleans
        self.periodicity = (False, False, False)
        #   self.current_time           <= simulation time in code units
        self.current_time = float(self._handle.time.values)
        #
        # We also set up cosmological information.  Set these to zero if
        # non-cosmological.
        #
        #   self.cosmological_simulation    <= int, 0 or 1
        #   self.current_redshift           <= float
        #   self.omega_lambda               <= float
        #   self.omega_matter               <= float
        #   self.hubble_constant            <= float
        self.cosmological_simulation = 0.0
        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        try:
            import xarray
        except ImportError:
            return False

        try:
            ds = xarray.open_dataset(args[0])
        except (FileNotFoundError, OSError):
            return False

        try:
            conventions = ds.attrs['Conventions']
        except KeyError:
            return False

        if 'CF/Radial' in conventions:
            return True

        return False
