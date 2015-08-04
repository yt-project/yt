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
import stat
import os
import time

from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.geometry.unstructured_mesh_handler import \
     UnstructuredIndex
from yt.data_objects.static_output import \
    Dataset
from yt.utilities.io_handler import \
    io_registry
from .io import \
    IOHandlerExodusII, mylog
from .fields import \
    ExodusIIFieldInfo
from .util import \
    load_info_records, sanitize_string

class ExodusIIGrid(AMRGridPatch):
    _id_offset = 1
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

class ExodusIIHierarchy(UnstructuredIndex):
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
        self.parameter_filename = filename
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
        self.dimensionality             = self.ds.variables['coor_names'].shape[0]
        self.parameters['info_records'] = self._load_info_records()
        self.unique_identifier          = self._get_unique_identifier()
        self.current_time               = self._get_current_time()
        self.parameters['num_elem']     = self.ds['eb_status'].shape[0]
        self.parameters['var_names']    = self._get_var_names()
        self.parameters['nod_names']    = self._get_nod_names()
        self.parameters['coordinates']  = self._load_coordinates()
        self.parameters['connectivity'] = self._load_connectivity()
        self.domain_left_edge           = self._load_domain_edge(0)
        self.domain_right_edge          = self._load_domain_edge(1)
        self.periodicity                = (False, False, False)
        self.cosmological_simulation    = 0
        self.current_redshift           = 0
        self.omega_lambda               = 0
        self.omega_matter               = 0
        self.hubble_constant            = 0

    def _load_info_records(self):
        try:
            return load_info_records(self.ds.variables['info_records'])
        except KeyError:
            mylog.warning("No info_records found")
            return []

    def _get_unique_identifier(self):
        try:
            return self.parameters['info_records']['Version Info']['Executable Timestamp']
        except KeyError:
            return self.parameter_filename.__hash__()

    def _get_current_time(self):
        try:
            return self.parameters['info_records']['Version Info']['Current Time']
        except KeyError:
            return 0.0

    def _get_var_names(self):
        try:
            return [sanitize_string(v.tostring()) for v in
                    self.ds.variables["name_elem_var"]]
        except KeyError:
            mylog.warning("name_elem_var not found")
            return []

    def _get_nod_names(self):
        try:
            return [sanitize_string(v.tostring()) for v in
                    self.ds.variables["name_nod_var"]]
        except KeyError:
            mylog.warning("name_nod_var not found")
            return []

    def _load_coordinates(self):
        if self.dimensionality == 3:
            coord_axes = 'xyz'
        elif self.dimensionality == 2:
            coord_axes = 'xy'

        mylog.info("Loading coordinates")

        if 'coord' in self.ds.variables.keys():
            return np.array([coord for coord in
                             self.ds.variables["coord"][:]]).transpose().copy()
        else:
            return np.array([self.ds.variables["coord%s" % ax][:]
                             for ax in coord_axes]).transpose().copy()
    
    def _load_connectivity(self):
        mylog.info("Loading connectivity")
        connectivity = []
        for i in range(self.parameters['num_elem']):
            connectivity.append(self.ds.variables["connect%d" % (i+1)][:].astype("i8"))
        return connectivity

    def _load_domain_edge(self, domain_idx):
        return np.array([self.parameters['coordinates'][:,domain_idx].min(),
                         self.parameters['coordinates'][:,domain_idx].max()],
                        'float64')

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        return False


