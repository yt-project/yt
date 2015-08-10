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
from yt.data_objects.unstructured_mesh import \
     UnstructuredMesh
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

class ExodusIIUnstructuredMesh(UnstructuredMesh):
    _index_offset = 1

    def __init__(self, *args, **kwargs):
        super(ExodusIIUnstructuredMesh, self).__init__(*args, **kwargs)

class ExodusIIUnstructuredIndex(UnstructuredIndex):
    def __init__(self, ds, dataset_type = 'exodus_ii'):
        super(ExodusIIUnstructuredIndex, self).__init__(ds, dataset_type)

class ExodusIIDataset(Dataset):
    _index_class = ExodusIIHierarchy
    _field_info_class = ExodusIIFieldInfo

    def __init__(self,
                 filename,
                 dataset_type='exodus_ii',
                 storage_filename=None,
                 units_override=None):

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
        self._load_variables()
        self.dimensionality             = self.parameters['coor_names'].shape[0]
        self.parameters['info_records'] = self._load_info_records()
        self.unique_identifier          = self._get_unique_identifier()
        self.current_time               = self._get_current_time()
        self.parameters['num_elem']     = self.parameters['eb_status'].shape[0]
        self.parameters['var_names']    = self._get_var_names()
        self.parameters['nod_names']    = self._get_nod_names()
        self.parameters['coordinates']  = self._load_coordinates()
        self.parameters['connectivity'] = self._load_connectivity()
        self.data                       = self._load_data()
        self.domain_left_edge           = self._load_domain_edge(0)
        self.domain_right_edge          = self._load_domain_edge(1)
        self.periodicity                = (False, False, False)
        self.cosmological_simulation    = 0
        self.current_redshift           = 0
        self.omega_lambda               = 0
        self.omega_matter               = 0
        self.hubble_constant            = 0

    def _load_variables(self):
        ds = IOHandlerExodusII(self.parameter_filename).ds
        for key in ds.variables.keys():
            self.parameters[key] = ds.variables[key]

    def _load_info_records(self):
        try:
            return load_info_records(self.parameters['info_records'])
        except (KeyError, TypeError):
            mylog.warning("No info_records found")
            return []

    def _get_unique_identifier(self):
        try:
            return self.parameters['info_records']['Version Info']['Executable Timestamp']
        except (KeyError, TypeError):
            return self.parameter_filename.__hash__()

    def _get_current_time(self):
        try:
            return self.parameters['info_records']['Version Info']['Current Time']
        except (KeyError, TypeError):
            return 0.0

    def _get_var_names(self):
        try:
            return [sanitize_string(v.tostring()) for v in
                    self.parameters["name_elem_var"]]
        except (KeyError, TypeError):
            mylog.warning("name_elem_var not found")
            return []

    def _get_nod_names(self):
        try:
            return [sanitize_string(v.tostring()) for v in
                    self.parameters["name_nod_var"]]
        except (KeyError, TypeError):
            mylog.warning("name_nod_var not found")
            return []

    def _load_coordinates(self):
        if self.dimensionality == 3:
            coord_axes = 'xyz'
        elif self.dimensionality == 2:
            coord_axes = 'xy'

        mylog.info("Loading coordinates")

        if 'coord' in self.parameters.keys():
            return np.array([coord for coord in
                             self.parameters["coord"][:]]).transpose().copy()
        else:
            return np.array([self.parameters["coord%s" % ax][:]
                             for ax in coord_axes]).transpose().copy()
    
    def _load_connectivity(self):
        mylog.info("Loading connectivity")
        connectivity = []
        for i in range(self.parameters['num_elem']):
            connectivity.append(self.parameters["connect%d" % (i+1)][:].astype("i8"))
        return connectivity

    def _load_data(self):
        data = []
        for i in range(self.parameters['num_elem']):
            ci = self.parameters['connectivity'][i]
            vals = {}

            for j, var_name in enumerate(self.parameters['var_names']):
                vals['gas', var_name] = self.parameters["vals_elem_var%seb%s" % (j+1, i+1)][:].astype("f8")[-1,:]

            for j, nod_name in enumerate(self.parameters['nod_names']):
                # We want just for this set of nodes all the node variables
                # Use (ci - 1) to get these values
                vals['gas', nod_name] = self.parameters["vals_nod_var%s" % (j+1)][:].astype("f8")[-1, ci - 1, ...]

            data.append(vals)

        return data

    def _load_domain_edge(self, domain_idx):
        return np.array([self.parameters['coordinates'][:,domain_idx].min(),
                         self.parameters['coordinates'][:,domain_idx].max()],
                        'float64')

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        return False


