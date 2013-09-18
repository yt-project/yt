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

import h5py
import os
import numpy as np
import weakref
from yt.funcs import mylog
from yt.data_objects.unstructured_mesh import \
           SemiStructuredMesh
from yt.geometry.unstructured_mesh_handler import \
           UnstructuredGeometryHandler
from yt.data_objects.static_output import \
           StaticOutput
from yt.utilities.io_handler import \
    io_registry
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion

from .fields import MoabFieldInfo, KnownMoabFields

def _get_convert(fname):
    def _conv(data):
        return data.convert(fname)
    return _conv

class MoabHex8Mesh(SemiStructuredMesh):
    _connectivity_length = 8
    _index_offset = 1

class MoabHex8Hierarchy(UnstructuredGeometryHandler):

    def __init__(self, pf, data_style='h5m'):
        self.parameter_file = weakref.proxy(pf)
        self.data_style = data_style
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self._fhandle = h5py.File(self.hierarchy_filename,'r')

        UnstructuredGeometryHandler.__init__(self, pf, data_style)

        self._fhandle.close()

    def _initialize_mesh(self):
        con = self._fhandle["/tstt/elements/Hex8/connectivity"][:]
        con = np.asarray(con, dtype="int64")
        coords = self._fhandle["/tstt/nodes/coordinates"][:]
        coords = np.asarray(coords, dtype="float64")
        self.meshes = [MoabHex8Mesh(0, self.hierarchy_filename, con,
                                    coords, self)]

    def _detect_fields(self):
        self.field_list = self._fhandle['/tstt/elements/Hex8/tags'].keys()

    def _count_grids(self):
        self.num_grids = 1 #self._fhandle['/grid_parent_id'].shape[0]

    def _setup_data_io(self):
        self.io = io_registry[self.data_style](self.parameter_file)

class MoabHex8StaticOutput(StaticOutput):
    _hierarchy_class = MoabHex8Hierarchy
    _fieldinfo_fallback = MoabFieldInfo
    _fieldinfo_known = KnownMoabFields
    periodicity = (False, False, False)

    def __init__(self, filename, data_style='moab_hex8',
                 storage_filename = None):
        StaticOutput.__init__(self, filename, data_style)
        self.storage_filename = storage_filename
        self.filename = filename
        self._handle = h5py.File(self.parameter_filename, "r")

    def _set_units(self):
        """Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['cm'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        for unit in mpc_conversion.keys():
            self.units[unit] = 1.0 * mpc_conversion[unit] / mpc_conversion["cm"]
        for unit in sec_conversion.keys():
            self.time_units[unit] = 1.0 / sec_conversion[unit]

    def _parse_parameter_file(self):
        self._handle = f = h5py.File(self.parameter_filename, "r")
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
        self.parameters['Time'] = 1.0 # Hardcode time conversion for now.
        self.parameters["HydroMethod"] = 0 # Hardcode for now until field staggering is supported.

    @classmethod
    def _is_valid(self, *args, **kwargs):
        fname = args[0]
        return fname.endswith('.h5m')

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]

class PyneMeshHex8Hierarchy(UnstructuredGeometryHandler):

    def __init__(self, pf, data_style='moab_hex8_pyne'):
        self.parameter_file = weakref.proxy(pf)
        self.data_style = data_style
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.getcwd()
        self.pyne_mesh = pf.pyne_mesh

        super(PyneMeshHex8Hierarchy, self).__init__(pf, data_style)

    def _initialize_mesh(self):
        from itaps import iBase, iMesh
        ent = self.pyne_mesh.structured_set.getEntities(iBase.Type.vertex)
        coords = self.pyne_mesh.mesh.getVtxCoords(ent).astype("float64")
        vind = self.pyne_mesh.structured_set.getAdjEntIndices(
            iBase.Type.region, iMesh.Topology.hexahedron,
            iBase.Type.vertex)[1].indices.data
        vind.shape = (vind.size/8.0, 8)
        con = vind.astype("int64")
        self.meshes = [MoabHex8Mesh(0, self.hierarchy_filename, con,
                                    coords, self)]

    def _detect_fields(self):
        self.field_list = self.pyne_mesh.mesh.getAllTags(
            self.pyne_mesh.mesh.rootSet)

    def _count_grids(self):
        self.num_grids = 1

    def _setup_data_io(self):
        self.io = io_registry[self.data_style](self.parameter_file)

class PyneMoabHex8StaticOutput(StaticOutput):
    _hierarchy_class = PyneMeshHex8Hierarchy
    _fieldinfo_fallback = MoabFieldInfo
    _fieldinfo_known = KnownMoabFields
    periodicity = (False, False, False)

    def __init__(self, pyne_mesh, data_style='moab_hex8_pyne',
                 storage_filename = None):
        filename = "pyne_mesh_" + str(id(pyne_mesh))
        self.pyne_mesh = pyne_mesh
        StaticOutput.__init__(self, str(filename), data_style)
        self.storage_filename = storage_filename
        self.filename = filename

    def _set_units(self):
        """Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['cm'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        for unit in mpc_conversion.keys():
            self.units[unit] = 1.0 * mpc_conversion[unit] / mpc_conversion["cm"]
        for unit in sec_conversion.keys():
            self.time_units[unit] = 1.0 / sec_conversion[unit]

    def _parse_parameter_file(self):
        from itaps import iBase
        ent = self.pyne_mesh.structured_set.getEntities(iBase.Type.vertex)
        coords = self.pyne_mesh.mesh.getVtxCoords(ent)
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
        self.parameters['Time'] = 1.0 # Hardcode time conversion for now.
        self.parameters["HydroMethod"] = 0 # Hardcode for now until field staggering is supported.

    @classmethod
    def _is_valid(self, *args, **kwargs):
        return False

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]

