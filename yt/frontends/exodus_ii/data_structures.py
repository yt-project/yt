"""
Exodus II data structures



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import numpy as np
import re

from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from .io import \
    IOHandlerExodusII
from .fields import \
    ExodusIIFieldInfo
from .util import load_info_records

class ExodusIIGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, id, index, level, start, dimensions):
        AMRGridPatch.__init__(self, id, filename=index.index_filename,
                              index=index)
        self.Parent = []
        self.Children = []
        self.Level = level
        self.start_index = start.copy()
        self.stop_index = self.start_index + dimensions
        self.ActiveDimensions = dimensions.copy()

    def __repr__(self):
        return "ExodusIIGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class ExodusIIHierarchy(GridIndex):
    grid = ExodusIIGrid

    def __init__(self, ds, dataset_type='exodus_ii'):
        self.dataset_type = dataset_type
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        # This needs to set a self.field_list that contains all the available,
        # on-disk fields.
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
        #
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

class ExodusIIDataset(Dataset):
    _index_class = ExodusIIHierarchy
    _field_info_class = ExodusIIFieldInfo

    def __init__(self,
                 filename,
                 dataset_type='exodus_ii',
                 storage_filename=None,
                 units_override=None):

        self.ds = IOHandlerExodusII(filename).ds
        self.fluid_types += ('exodus_ii',)
        self.parameter_filename = storage_filename
        Dataset.__init__(self, filename, dataset_type,
                         units_override=units_override)
        self.storage_filename = storage_filename

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
        self.parameters['info_records'] = load_info_records(self.ds.variables['info_records'])
        self.unique_identifier          = self.parameters['info_records']['Version Info']['Executable Timestamp']
        self.current_time               = self.parameters['info_records']['Version Info']['Current Time']
        self.dimensionality             = self.ds.variables['coor_names'].shape[0]
        # self.domain_left_edge           = np.array([self.ds.coordinates[:,0].min(),
        #                                             self.ds.coordinates[:,0].max()],
        #                                            'float64')
        # self.domain_right_edge           = np.array([self.ds.coordinates[:,1].min(),
        #                                              self.ds.coordinates[:,1].max()],
        #                                             'float64')
        self.periodicity                = (False, False, False)
        self.cosmological_simulation    = 0
        self.current_redshift           = 0
        self.omega_lambda               = 0
        self.omega_matter               = 0
        self.hubble_constant            = 0

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        return False


