"""
Denovo data structures



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import numpy as np
import weakref
from yt.utilities.on_demand_imports import _h5py as h5py
from yt.utilities.logger import ytLogger as mylog

from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset

from .fields import DenovoFieldInfo


class DenovoGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        AMRGridPatch.__init__(self, id, filename=index.index_filename,
                              index=index)
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "DenovoGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class DenovoHierarchy(GridIndex):
    grid = DenovoGrid

    def __init__(self, ds, dataset_type='denovo'):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        GridIndex.__init__(self, ds, dataset_type)

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


class DenovoDataset(Dataset):
    _index_class = DenovoHierarchy
    _field_info_class = DenovoFieldInfo
    _suffix = ".out.h5"

    def __init__(self, filename,
                 dataset_type='denovo',
                 storage_filename=None,
                 units_override=None):

        self.fluid_types += ('denovo',)
        self._handle = HDF5FileHandler(filename)
        self.dataset_type = dataset_type

        self.geometry = 'cartesian'

        super(DenovoDataset, self).__init__(self, filename, dataset_type,
                         units_override=units_override)
        self.storage_filename = storage_filename
        self.cosmological_simulation = False

        self.parameters["HydroMethod"] = 'denovo'
        self.parameters
        # refinement factor between a grid and its subgrid
        # self.refine_by = 2

    def _set_code_unit_attributes(self):
        #
        # For now set the length mass and time units to what we expect.
        # Denovo does not currently output what units the flux is in.
        #
        #
        setdefaultattr(self, 'length_unit', self.quan(length_unit, "cm"))
        setdefaultattr(self, 'mass_unit', self.quan(length_unit, "g"))
        setdefaultattr(self, 'time_unit', self.quan(length_unit, "s"))
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "cm/s")
        # self.magnetic_unit = self.quan(1.0, "gauss")
        pass

    def _parse_parameter_file(self):
        #
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        self.parameters             <= full of code-specific items of use
        self.domain_left_edge       <= array of float64
        self.domain_right_edge      <= array of float64
        self.dimensionality         <= int
        self.domain_dimensions      <= array of int64

        # We know that Denovo datasets are never periodic, so periodicity will
        # be set to false.
        self.periodicity = (False, False, False)

        # There is no time-dependence in the denovo solutions at this time, so
        # this will be set to 0.0
        self.current_time = 0.0
        #
        # The next few paramaters are set to 0 because Denovo is a
        # non-cosmological simulation tool.
        #
        self.cosmological_simulation = 0
        self.current_redshift = 0
        self.omega_lambda = 0
        self.omega_matter = 0
        self.hubble_constant = 0
        pass

    @classmethod
    def _is_valid(self, *args, **kwargs):

        warn_h5py(args[0])

        try:
            fileh = h5py.File(args[0], 'r')

            # for now check that denovo is in the solution file. This will need
            # to be updated for fwd/adjoint runs in the future with multiple
            # arguments for a valid dataset.
            valid = "denovo" in fileh["/"]
            fileh.close()
            return valid
        except:
            pass
        return False
