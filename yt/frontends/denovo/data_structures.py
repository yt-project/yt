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
from yt.data_objects.unstructured_mesh import \
    SemiStructuredMesh
from yt.geometry.unstructured_mesh_handler import \
    UnstructuredIndex
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset

from .fields import DenovoFieldInfo


# class DenovoGrid(AMRGridPatch):
#     _id_offset = 0
#
#     def __init__(self, id, index, level):
#         super(DenovoGrid, self).__init__(id, filename=index.index_filename,
#                               index=index)
#         self.Parent = None
#         self.Children = []
#         self.Level = level
#
#     def __repr__(self):
#         return "DenovoGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class DenovoMesh(SemiStructuredMesh):
    _connectivity_length = 8
    _index_offset = 1

class DenovoHierarchy(UnstructuredIndex):

    def __init__(self, ds, dataset_type='denovo'):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)

        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self._fhandle = h5py.File(self.index_filename, 'r')

        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super(DenovoHierarchy, self).__init__(ds, dataset_type)

    def _initialize_mesh(self):
        from yt.frontends.stream.data_structures import hexahedral_connectivity
        coords, conn = hexahedral_connectivity(ds.paramaters.mesh_x,
                ds.parameters.mesh_y, ds.parameters.mesh_z)
        self.meshes = [DenovoMesh(0, self.index_filename, conn, coords, self)]


    def _detect_output_fields(self):

        self.field_list = [("denovo", key) for key, val in
                self._fhandle['/denovo/db/hdf5_db'].items() if val.value]

        # right now we don't want to read in a few of these fields, so we'll
        # manually delete them here until later.
        #
        if ('denovo','angular_mesh') in self.field_list:
            self.field_list.remove('denovo','angular_mesh')
        self.field_list.remove('denovo','block')

    def _count_grids(self):
        self.num_grids=1

    # not sure if this is needed since it's not in the hexahedral OR Moab
    # Hierarchy classes.
    # def _parse_index(self):
    #     # This needs to fill the following arrays, where N is self.num_grids:
    #     #   self.grid_left_edge         (N, 3) <= float64
    #     #   self.grid_right_edge        (N, 3) <= float64
    #     #   self.grid_dimensions        (N, 3) <= int
    #     #   self.grid_particle_count    (N, 1) <= int
    #     #   self.grid_levels            (N, 1) <= int
    #     #   self.grids                  (N, 1) <= grid objects
    #     #   self.max_level = self.grid_levels.max()
    #     pass

    # def _populate_grid_objects(self):
    #     # For each grid, this must call:
    #     #   grid._prepare_grid()
    #     #   grid._setup_dx()
    #     # This must also set:
    #     #   grid.Children <= list of child grids
    #     #   grid.Parent   <= parent grid
    #     # This is handled by the frontend because often the children must be
    #     # identified.
    #     pass


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

        super(DenovoDataset, self).__init__(filename, dataset_type,
                         units_override=units_override)
        self.storage_filename = storage_filename
        self.filename=filename

        # Note for later: if this is set in _parse_parameter_file, does it need
        # to be set here?
        self.cosmological_simulation = False

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
        pass

    def _parse_parameter_file(self):
        #
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        self._handle = h5py.File(self.parameter_filename, "r")
        self.parameters = self._load_parameters()
        self.domain_left_edge, self.domain_right_edge = self._load_domain_edge()
        self.domain_dimensions = self.domain_right_edge - self.domain_left_edge
        self.dimensionality = len(self.domain_dimensions)

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

    def _load_parameters():
        """

        loads in all of the parameters located in the 'denovo' group that are
        not groups themselves and are not flux or source fields.

        """

        for key,val in self._handle["/denovo"].items():
            if isinstance(val, h5py.Dataset):
                # skip the fields that will be stored elsewhere
                if any([key == 'flux', key == 'source']):
                    pass
                else:
                    self.parameters[key]=val.value

        return self.parameters


    def _load_domain_edge():
        """

        Loads the boundaries for the domain edge using the min and max values
        obtained from the x, y, and z meshes.

        """

        mylog.info("calculating domain boundaries")

        if 'mesh_x' in self.parameters:
            left_edge = [self.parameters['mesh_x'].min(),
            self.parameters['mesh_y'].min(),
            self.parameters['mesh_z'].min()]

            right_edge = [self.parameters['mesh_x'].max(),
            self.parameters['mesh_y'].max(),
            self.parameters['mesh_z'].max()]

        else:
            left_edge = []
            right_edge = []

        return left_edge, right_edge


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
