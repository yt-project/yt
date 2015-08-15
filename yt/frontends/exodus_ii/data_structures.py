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

from yt.geometry.unstructured_mesh_handler import \
    UnstructuredIndex
from yt.data_objects.unstructured_mesh import \
    UnstructuredMesh
from yt.data_objects.static_output import \
    Dataset
from yt.utilities.io_handler import \
    io_registry
from .io import \
    IOHandlerExodusII, \
    NetCDF4FileHandler
from yt.utilities.logger import ytLogger as mylog
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
        self._connectivity_length = self.meshes[0].connectivity_indices.shape[1]

    def _initialize_mesh(self):
        self.meshes = [ExodusIIUnstructuredMesh(
            mesh_id, self.index_filename, conn_ind, conn_coord, self)
                       for mesh_id, (conn_ind, conn_coord) in
                       enumerate(zip(self.dataset.parameters['connectivity'],
                                     self.dataset.parameters['coordinates']))]

    def _setup_data_io(self):
        self.io = io_registry[self.dataset_type](self.ds)

    def _detect_output_fields(self):
        ftype = 'exodus_ii'
        fnames = self.dataset.parameters['nod_names']
        self.field_list = [(ftype, fname) for fname in fnames]


class ExodusIIDataset(Dataset):
    _index_class = ExodusIIUnstructuredIndex
    _field_info_class = ExodusIIFieldInfo

    def __init__(self,
                 filename,
                 dataset_type='exodus_ii',
                 storage_filename=None,
                 units_override=None):

        self.fluid_types               += ('exodus_ii',)
        self.parameter_filename         = filename
        Dataset.__init__(self, filename, dataset_type,
                         units_override = units_override)
        self.index_filename             = filename
        self.storage_filename           = storage_filename

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        self.length_unit = self.quan(1.0, "cm")
        self.mass_unit = self.quan(1.0, "g")
        self.time_unit = self.quan(1.0, "s")
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "cm/s")
        # self.magnetic_unit = self.quan(1.0, "gauss")

    def _parse_parameter_file(self):
        self._load_variables()
        self._read_glo_var()
        self.dimensionality             = self.parameters['coor_names'].shape[0]
        self.parameters['info_records'] = self._load_info_records()
        self.unique_identifier          = self._get_unique_identifier()
        self.current_time               = self._get_current_time()
        self.parameters['num_meshes']   = self.parameters['eb_status'].shape[0]
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
        self.refine_by                  = 0

    def _load_variables(self):
        """
        Loads each key-pair in the Exodus II input file
        as a parameter
        """
        handle = NetCDF4FileHandler(self.parameter_filename).dataset
        for key in handle.variables.keys():
            self.parameters[key] = handle.variables[key]

    def _read_glo_var(self):
        """
        Adds each global variable to the dict of parameters

        """
        names = self._get_glo_names()
        values = self.parameters['vals_glo_var'][:].transpose() 
        for name, value in zip(names, values):
            self.parameters[name] = value

    def _load_info_records(self):
        """
        Returns parsed version of the info_records.
        """
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

    def _get_glo_names(self):
        """
        Returns the names of the global vars, if available
        """
        try:
            return [sanitize_string(v.tostring()) for v in
                    self.parameters["name_glo_var"]]
        except (KeyError, TypeError):
            mylog.warning("name_glo_var not found")
            return []

    def _get_var_names(self):
        """
        Returns the names of the element vars, if available
        """
        try:
            return [sanitize_string(v.tostring()) for v in
                    self.parameters["name_elem_var"]]
        except (KeyError, TypeError):
            mylog.warning("name_elem_var not found")
            return []

    def _get_nod_names(self):
        """
        Returns the names of the node vars, if available
        """
        try:
            return [sanitize_string(v.tostring()) for v in
                    self.parameters["name_nod_var"]]
        except (KeyError, TypeError):
            mylog.warning("name_nod_var not found")
            return []

    def _load_coordinates(self):
        """
        Loads the coordinates for the mesh
        """
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
        """
        Loads the connectivity data for the mesh
        """
        mylog.info("Loading connectivity")
        connectivity = []
        for i in range(self.parameters['num_meshes']):
            connectivity.append(self.parameters["connect%d" % (i+1)][:].astype("i8"))
        return connectivity

    def _load_data(self):
        """
        Loads the fluid data
        """
        data = []
        for i in range(self.parameters['num_meshes']):
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
        """
        Loads the boundaries for the domain edge

        Parameters:
        - domain_idx: 0 corresponds to the left edge, 1 corresponds to the right edge
        """
        if domain_idx == 0:
            return self.parameters['coordinates'].min(axis=0)
        if domain_idx == 1:
            return self.parameters['coordinates'].max(axis=0)        

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            filename = args[0]
            fhandle = NetCDF4FileHandler(filename)
            fmt = fhandle.dataset.file_format
            return fmt.upper().startswith('NETCDF')
        except:
            pass
