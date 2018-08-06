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

import os
import numpy as np
import weakref

from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from .fields import CFRadialFieldInfo


class CFRadialGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        super(CFRadialGrid, self).__init__(
            self, id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level

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
        super(CFRadialHierarchy, self).__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        # This needs to set a self.field_list that contains all the available,
        # on-disk fields. No derived fields should be defined here.
        # NOTE: Each should be a tuple, where the first element is the on-disk
        # fluid type or particle type.  Convention suggests that the on-disk
        # fluid type is usually the dataset_type and the on-disk particle type
        # (for a single population of particles) is "io".
        pass

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
        super(CFRadialDataset, self).__init__(self, filename, dataset_type,
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
        # self.mass_unit = self.quan(1.0, "g")
        # self.time_unit = self.quan(1.0, "s")
        # self.time_unit = self.quan(1.0, "s")
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "cm/s")
        # self.magnetic_unit = self.quan(1.0, "gauss")
        pass

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
        pass

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
        except FileNotFoundError:
            return False

        try:
            conventions = ds.attrs['Conventions']
        except KeyError:
            return False

        if 'CF/Radial' in conventions:
            return True

        return False
