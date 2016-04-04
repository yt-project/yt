"""
Data structures for Athena.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import os
import weakref
import glob
import lxml

from yt.funcs import \
    mylog, \
    ensure_tuple
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from yt.utilities.lib.misc_utilities import \
    get_box_grids_level
from yt.geometry.geometry_handler import \
    YTDataChunk
from yt.utilities.file_handler import \
    HDF5FileHandler

from .fields import AthenaPPFieldInfo

geometry_map = {"3DRectMesh": "cartesian"}

class AthenaPPGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        AMRGridPatch.__init__(self, id, filename = index.index_filename,
                              index = index)
        self.Parent = []
        self.Children = []
        self.Level = level

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        if len(self.Parent) > 0:
            self.dds = self.Parent[0].dds / self.ds.refine_by
        else:
            LE, RE = self.index.grid_left_edge[id,:], \
                     self.index.grid_right_edge[id,:]
            self.dds = self.ds.arr((RE-LE)/self.ActiveDimensions, "code_length")
        if self.ds.dimensionality < 2: self.dds[1] = 1.0
        if self.ds.dimensionality < 3: self.dds[2] = 1.0
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

    def __repr__(self):
        return "AthenaPPGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class AthenaPPHierarchy(GridIndex):

    grid = AthenaPPGrid
    _dataset_type='athena++'
    _data_file = None

    def __init__(self, ds, dataset_type='athena++'):
        self.dataset = weakref.proxy(ds)
        self.directory = os.path.dirname(self.dataset.filename)
        self.dataset_type = dataset_type
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.datafile
        self._handle = ds._handle
        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        self._field_map = dict((("athena++", k), v) for k, v in self.dataset._field_map.items())
        self.field_list = list(self._field_map.keys())

    def _count_grids(self):
        self.num_grids = self._handle.attrs["TotalMeshBlock"][0]

    def _parse_index(self):
        coord_fields = self.dataset._coord_fields
        num_grids = self._handle.attrs["TotalMeshBlock"][0]

        self.grid_left_edge = np.zeros((num_grids, 3), dtype='float64')
        self.grid_right_edge = np.zeros((num_grids, 3), dtype='float64')
        self.grid_dimensions = np.zeros((num_grids, 3), dtype='int32')
        dds = np.zeros((num_grids, 3), dtype='float64')

        for i in range(num_grids):
            x = self._handle["MeshBlock%d" % i][coord_fields[0]]
            y = self._handle["MeshBlock%d" % i][coord_fields[1]]
            z = self._handle["MeshBlock%d" % i][coord_fields[2]]
            dds[i] = np.array([x[1]-x[0], y[1]-y[0], z[1]-z[0]], dtype='float64')
            self.grid_left_edge[i] = np.array([x[0], y[0], z[0]], dtype='float64')
            self.grid_dimensions[i] = self._handle.attrs["MeshBlockSize"]
        self.grid_right_edge = self.grid_left_edge + self.grid_dimensions*dds

        new_dle = np.min(self.grid_left_edge, axis=0)
        new_dre = np.max(self.grid_right_edge, axis=0)
        self.dataset.domain_left_edge[:] = np.round(new_dle, decimals=12)[:]
        self.dataset.domain_right_edge[:] = np.round(new_dre, decimals=12)[:]
        self.dataset.domain_width = \
                (self.dataset.domain_right_edge -
                 self.dataset.domain_left_edge)
        self.dataset.domain_center = \
                0.5*(self.dataset.domain_left_edge +
                     self.dataset.domain_right_edge)

        dx_root = (self.dataset.domain_right_edge-
                   self.dataset.domain_left_edge)/self.dataset.domain_dimensions

        # This next line assumes refine_by == 2!!
        levels = np.round(np.log2(dds[:,0]/dx_root[0])).astype('int')

        self.grid_left_edge = self.ds.arr(self.grid_left_edge, "code_length")
        self.grid_right_edge = self.ds.arr(self.grid_right_edge, "code_length")

        self.grids = np.empty(self.num_grids, dtype='object')
        for i in range(num_grids):
            self.grids[i] = self.grid(i, self, levels[i])

        if self.dataset.dimensionality <= 2:
            self.grid_right_edge[:,2] = self.dataset.domain_right_edge[2]
        if self.dataset.dimensionality == 1:
            self.grid_right_edge[:,1:] = self.dataset.domain_right_edge[1:]
        self.grid_particle_count = np.zeros([self.num_grids, 1], dtype='int64')

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()
        self._reconstruct_parent_child()
        self.max_level = self._handle.attrs["MaxLevel"][0]

    def _reconstruct_parent_child(self):
        mask = np.empty(len(self.grids), dtype='int32')
        mylog.debug("First pass; identifying child grids")
        for i, grid in enumerate(self.grids):
            get_box_grids_level(self.grid_left_edge[i,:],
                                self.grid_right_edge[i,:],
                                self.grid_levels[i] + 1,
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask)
            grid.Children = [g for g in self.grids[mask.astype("bool")] if g.Level == grid.Level + 1]
        mylog.debug("Second pass; identifying parents")
        for i, grid in enumerate(self.grids): # Second pass
            for child in grid.Children:
                child.Parent.append(grid)

    def _get_grid_children(self, grid):
        mask = np.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        return [g for g in self.grids[mask] if g.Level == grid.Level + 1]

    def _chunk_io(self, dobj, cache = True, local_only = False):
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in gobjs:
            yield YTDataChunk(dobj, "io", [subset],
                              self._count_selection(dobj, [subset]),
                              cache = cache)

class AthenaPPDataset(Dataset):
    _index_class = AthenaPPHierarchy
    _field_info_class = AthenaPPFieldInfo
    _dataset_type = "athena++"

    def __init__(self, filename, dataset_type='athena++',
                 storage_filename=None, parameters=None,
                 units_override=None, unit_system="cgs"):
        self.fluid_types += ("athena++",)
        if self._handle is not None: return
        if parameters is None:
            parameters = {}
        self.specified_parameters = parameters
        if units_override is None:
            units_override = {}
        Dataset.__init__(self, filename, dataset_type, units_override=units_override,
                         unit_system=unit_system)
        self.filename = filename
        if storage_filename is None:
            storage_filename = '%s.yt' % filename.split('/')[-1]
        self.storage_filename = storage_filename
        self.backup_filename = self.filename[:-4] + "_backup.gdf"
        # Unfortunately we now have to mandate that the index gets
        # instantiated so that we can make sure we have the correct left
        # and right domain edges.
        self.index

    def _set_code_unit_attributes(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        if "length_unit" not in self.units_override:
            self.no_cgs_equiv_length = True
        for unit, cgs in [("length", "cm"), ("time", "s"), ("mass", "g")]:
            # We set these to cgs for now, but they may be overridden later.
            mylog.warning("Assuming 1.0 = 1.0 %s", cgs)
            setattr(self, "%s_unit" % unit, self.quan(1.0, cgs))

    def set_code_units(self):
        super(AthenaPPDataset, self).set_code_units()
        mag_unit = getattr(self, "magnetic_unit", None)
        vel_unit = getattr(self, "velocity_unit", None)
        if mag_unit is None:
            self.magnetic_unit = np.sqrt(4*np.pi * self.mass_unit /
                                         (self.time_unit**2 * self.length_unit))
        self.magnetic_unit.convert_to_units("gauss")
        self.unit_registry.modify("code_magnetic", self.magnetic_unit)
        if vel_unit is None:
            self.velocity_unit = self.length_unit/self.time_unit
        self.unit_registry.modify("code_velocity", self.velocity_unit)

    def _parse_parameter_file(self):
        tree = lxml.etree.parse(self.parameter_filename)
        root = tree.getroot()
        grids = root.findall("./Domain/Grid")[0]
        grid0 = grids[0]
        geometry = grid0.find("Topology").attrib["TopologyType"]
        coords = grid0.find("Geometry")
        self._coord_fields = [coord.text.split("/")[-1] for coord in coords]
        self._field_map = {}
        for field in grid0.findall("Attribute"):
            self._field_map[field.attrib["Name"]] = field[0].text.split("/")[-1]
        self.datafile = coords[0].text.split(":")[0]
        self._handle = HDF5FileHandler(self.datafile)

        xmin = self._handle["MeshBlock0"][self._coord_fields[0]][0]
        ymin = self._handle["MeshBlock0"][self._coord_fields[1]][0]
        zmin = self._handle["MeshBlock0"][self._coord_fields[2]][0]
        mylog.info("Temporarily setting a likely bogus domain_right_edge and "
                   "domain_left_edge. This will be corrected later.")
        self.domain_left_edge = np.array([xmin, ymin, zmin], dtype='float64')
        self.domain_right_edge = -self.domain_left_edge
        self.domain_width = self.domain_right_edge-self.domain_left_edge
        self.domain_dimensions = self._handle.attrs["RootGridSize"]

        self.refine_by = 2
        dimensionality = 3
        if self.domain_dimensions[2] == 1 :
            dimensionality = 2
        if self.domain_dimensions[1] == 1 :
            dimensionality = 1
        self.dimensionality = dimensionality
        self.current_time = self._handle.attrs["Time"][0]
        self.unique_identifier = self.parameter_filename.__hash__()
        self.cosmological_simulation = False
        self.num_ghost_zones = 0
        self.field_ordering = 'fortran'
        self.boundary_conditions = [1]*6
        if 'periodicity' in self.specified_parameters:
            self.periodicity = ensure_tuple(self.specified_parameters['periodicity'])
        else:
            self.periodicity = (True,True,True,)
        if 'gamma' in self.specified_parameters:
            self.gamma = float(self.specified_parameters['gamma'])
        else:
            self.gamma = 5./3.

        self.current_redshift = self.omega_lambda = self.omega_matter = \
            self.hubble_constant = self.cosmological_simulation = 0.0
        self.parameters['Time'] = self.current_time # Hardcode time conversion for now.
        self.parameters["HydroMethod"] = 0 # Hardcode for now until field staggering is supported.
        if "gamma" in self.specified_parameters:
            self.parameters["Gamma"] = self.specified_parameters["gamma"]
        else:
            self.parameters["Gamma"] = 5./3.
        self.geometry = self.geometry_map.get(geometry, "cartesian")

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            if 'athdf5.xdmf' in args[0]:
                return True
        except:
            pass
        return False

    @property
    def _skip_cache(self):
        return True

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]
