"""
Data structures for MOAB Hex8.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import os
import numpy as np
import weakref
from yt.data_objects.unstructured_mesh import \
    SemiStructuredMesh
from yt.funcs import \
    setdefaultattr
from yt.geometry.unstructured_mesh_handler import \
    UnstructuredIndex
from yt.data_objects.static_output import \
    Dataset
from yt.utilities.file_handler import HDF5FileHandler

from .fields import MoabFieldInfo, PyneFieldInfo

class MoabHex8Mesh(SemiStructuredMesh):
    _connectivity_length = 8
    _index_offset = 1

class MoabHex8Hierarchy(UnstructuredIndex):

    def __init__(self, ds, dataset_type='h5m'):
        self.dataset = weakref.proxy(ds)
        self.dataset_type = dataset_type
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self._fhandle = h5py.File(self.index_filename,'r')

        UnstructuredIndex.__init__(self, ds, dataset_type)

        self._fhandle.close()

    def _initialize_mesh(self):
        con = self._fhandle["/tstt/elements/Hex8/connectivity"][:]
        con = np.asarray(con, dtype="int64")
        coords = self._fhandle["/tstt/nodes/coordinates"][:]
        coords = np.asarray(coords, dtype="float64")
        self.meshes = [MoabHex8Mesh(0, self.index_filename, con,
                                    coords, self)]

    def _detect_output_fields(self):
        self.field_list = [("moab", f) for f in 
            self._fhandle['/tstt/elements/Hex8/tags'].keys()]

    def _count_grids(self):
        self.num_grids = 1

class MoabHex8Dataset(Dataset):
    _index_class = MoabHex8Hierarchy
    _field_info_class = MoabFieldInfo
    periodicity = (False, False, False)

    def __init__(self, filename, dataset_type='moab_hex8',
                 storage_filename = None, units_override=None,
                 unit_system="cgs"):
        self.fluid_types += ("moab",)
        Dataset.__init__(self, filename, dataset_type,
                         units_override=units_override,
                         unit_system=unit_system)
        self.storage_filename = storage_filename
        self.filename = filename
        self._handle = HDF5FileHandler(filename)

    def _set_code_unit_attributes(self):
        # Almost everything is regarded as dimensionless in MOAB, so these will
        # not be used very much or at all.
        setdefaultattr(self, 'length_unit', self.quan(1.0, "cm"))
        setdefaultattr(self, 'time_unit', self.quan(1.0, "s"))
        setdefaultattr(self, 'mass_unit', self.quan(1.0, "g"))

    def _parse_parameter_file(self):
        self._handle = h5py.File(self.parameter_filename, "r")
        coords = self._handle["/tstt/nodes/coordinates"]
        self.domain_left_edge = coords[0]
        self.domain_right_edge = coords[-1]
        self.domain_dimensions = self.domain_right_edge - self.domain_left_edge
        self.refine_by = 2
        self.dimensionality = len(self.domain_dimensions)
        self.current_time = 0.0
        self.unique_identifier = self.parameter_filename
        self.cosmological_simulation = False
        self.num_ghost_zones = 0
        self.current_redshift = self.omega_lambda = self.omega_matter \
                              = self.hubble_constant \
                              = self.cosmological_simulation = 0.0

    @classmethod
    def _is_valid(self, *args, **kwargs):
        fname = args[0]
        return fname.endswith('.h5m')

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]

class PyneHex8Mesh(SemiStructuredMesh):
    _connectivity_length = 8
    _index_offset = 0

class PyneMeshHex8Hierarchy(UnstructuredIndex):

    def __init__(self, ds, dataset_type='moab_hex8_pyne'):
        self.dataset = weakref.proxy(ds)
        self.dataset_type = dataset_type
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.getcwd()
        self.pyne_mesh = ds.pyne_mesh

        super(PyneMeshHex8Hierarchy, self).__init__(ds, dataset_type)

    def _initialize_mesh(self):
        from itaps import iBase, iMesh
        
        ents = list(self.pyne_mesh.structured_iterate_vertex())
        coords = self.pyne_mesh.mesh.getVtxCoords(ents).astype("float64")
        vind = self.pyne_mesh.mesh.rootSet.getAdjEntIndices(
            iBase.Type.region, iMesh.Topology.hexahedron,
            iBase.Type.vertex)[1].indices.data.astype("int64")
        # Divide by float so it throws an error if it's not 8
        vind.shape = (vind.shape[0] / 8.0, 8)
        self.meshes = [PyneHex8Mesh(0, self.index_filename,
                                    vind, coords, self)]

    def _detect_output_fields(self):
        self.field_list = [("pyne", f) for f in self.pyne_mesh.tags.keys()]

    def _count_grids(self):
        self.num_grids = 1

class PyneMoabHex8Dataset(Dataset):
    _index_class = PyneMeshHex8Hierarchy
    _fieldinfo_fallback = MoabFieldInfo
    _field_info_class = PyneFieldInfo
    periodicity = (False, False, False)

    def __init__(self, pyne_mesh, dataset_type='moab_hex8_pyne',
                 storage_filename = None, units_override=None,
                 unit_system="cgs"):
        self.fluid_types += ("pyne",)
        filename = "pyne_mesh_" + str(id(pyne_mesh))
        self.pyne_mesh = pyne_mesh
        Dataset.__init__(self, str(filename), dataset_type,
                         units_override=units_override,
                         unit_system=unit_system)
        self.storage_filename = storage_filename
        self.filename = filename

    def _set_code_unit_attributes(self):
        # Almost everything is regarded as dimensionless in MOAB, so these will
        # not be used very much or at all.
        setdefaultattr(self, 'length_unit', self.quan(1.0, "cm"))
        setdefaultattr(self, 'time_unit', self.quan(1.0, "s"))
        setdefaultattr(self, 'mass_unit', self.quan(1.0, "g"))

    def _parse_parameter_file(self):
        #  not sure if this import has side-effects so I'm not deleting it
        from itaps import iBase  # NOQA

        ents = list(self.pyne_mesh.structured_iterate_vertex())
        coords = self.pyne_mesh.mesh.getVtxCoords(ents)
        self.domain_left_edge = coords[0]
        self.domain_right_edge = coords[-1]
        self.domain_dimensions = self.domain_right_edge - self.domain_left_edge
        self.refine_by = 2
        self.dimensionality = len(self.domain_dimensions)
        self.current_time = 0.0
        self.unique_identifier = self.parameter_filename
        self.cosmological_simulation = False
        self.num_ghost_zones = 0
        self.current_redshift = self.omega_lambda = self.omega_matter \
                              = self.hubble_constant \
                              = self.cosmological_simulation = 0.0

    @classmethod
    def _is_valid(self, *args, **kwargs):
        return False

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]

