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

from yt.funcs import \
    mylog, get_pbar, \
    ensure_tuple
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from yt.geometry.geometry_handler import \
    YTDataChunk
from yt.utilities.file_handler import \
    HDF5FileHandler
from yt.geometry.unstructured_mesh_handler import \
    UnstructuredIndex
from yt.data_objects.unstructured_mesh import \
    SemiStructuredMesh
from itertools import chain, product
from .fields import AthenaPPFieldInfo

geom_map = {"cartesian": "cartesian",
            "cylindrical": "cylindrical",
            "spherical_polar": "spherical",
            "minkowski": "cartesian",
            "tilted": "cartesian",
            "sinusoidal": "cartesian",
            "schwarzschild": "spherical",
            "kerr-schild": "spherical"}

_cis = np.fromiter(chain.from_iterable(product([0,1], [0,1], [0,1])),
                   dtype=np.int64, count = 8*3)
_cis.shape = (8, 3)

class AthenaPPLogarithmicMesh(SemiStructuredMesh):
    _index_offset = 0

    def __init__(self, mesh_id, filename, connectivity_indices,
                 connectivity_coords, index, blocks, dims):
        super(AthenaPPLogarithmicMesh, self).__init__(mesh_id, filename, 
                                                      connectivity_indices,
                                                      connectivity_coords, index)
        self.mesh_blocks = blocks
        self.mesh_dims = dims

class AthenaPPLogarithmicIndex(UnstructuredIndex):
    def __init__(self, ds, dataset_type = 'athena_pp'):
        self._handle = ds._handle
        super(AthenaPPLogarithmicIndex, self).__init__(ds, dataset_type)
        self.index_filename = self.dataset.filename
        self.directory = os.path.dirname(self.dataset.filename)
        self.dataset_type = dataset_type

    def _initialize_mesh(self):
        mylog.debug("Setting up meshes.")
        num_blocks = self._handle.attrs["NumMeshBlocks"]
        log_loc = self._handle['LogicalLocations']
        levels = self._handle["Levels"]
        x1f = self._handle["x1f"]
        x2f = self._handle["x2f"]
        x3f = self._handle["x3f"]
        nbx, nby, nbz = tuple(np.max(log_loc, axis=0)+1)
        nlevel = self._handle.attrs["MaxLevel"]+1

        nb = np.array([nbx, nby, nbz], dtype='int')
        self.mesh_factors = np.ones(3, dtype='int')*((nb > 1).astype("int")+1)

        block_grid = -np.ones((nbx,nby,nbz,nlevel), dtype=np.int)
        block_grid[log_loc[:,0],log_loc[:,1],log_loc[:,2],levels[:]] = np.arange(num_blocks)

        block_list = np.arange(num_blocks, dtype='int64')
        bc = []
        for i in range(num_blocks):
            if block_list[i] >= 0:
                ii, jj, kk = log_loc[i]
                neigh = block_grid[ii:ii+2,jj:jj+2,kk:kk+2,levels[i]]
                if np.all(neigh > -1):
                    loc_ids = neigh.transpose().flatten()
                    bc.append(loc_ids)
                    block_list[loc_ids] = -1
                else:
                    bc.append(np.array(i))
                    block_list[i] = -1

        num_meshes = len(bc)

        self.meshes = []
        pbar = get_pbar("Constructing meshes", num_meshes)
        for i in range(num_meshes):
            ob = bc[i][0]
            x = x1f[ob,:]
            y = x2f[ob,:]
            z = x3f[ob,:]
            if nbx > 1:
                x = np.concatenate([x, x1f[bc[i][1],1:]])
            if nby > 1:
                y = np.concatenate([y, x2f[bc[i][2],1:]])
            if nbz > 1:
                z = np.concatenate([z, x3f[bc[i][4],1:]])
            nxm = x.size
            nym = y.size
            nzm = z.size
            coords = np.zeros((nxm, nym, nzm, 3), dtype="float64", order="C")
            coords[:,:,:,0] = x[:,None,None]
            coords[:,:,:,1] = y[None,:,None]
            coords[:,:,:,2] = z[None,None,:]
            coords.shape = (nxm * nym * nzm, 3)
            cycle = np.rollaxis(np.indices((nxm-1,nym-1,nzm-1)), 0, 4)
            cycle.shape = ((nxm-1)*(nym-1)*(nzm-1), 3)
            off = _cis + cycle[:, np.newaxis]
            connectivity = ((off[:,:,0] * nym) + off[:,:,1]) * nzm + off[:,:,2]
            mesh = AthenaPPLogarithmicMesh(i, self.index_filename, connectivity,
                                           coords, self, bc[i],
                                           np.array([nxm-1, nym-1, nzm-1]))
            self.meshes.append(mesh)
            pbar.update(i)
        pbar.finish()
        mylog.debug("Done setting up meshes.")

    def _detect_output_fields(self):
        self.field_list = [("athena_pp", k) for k in self.ds._field_map]

    def _chunk_io(self, dobj, cache = True, local_only = False):
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in gobjs:
            yield YTDataChunk(dobj, "io", [subset],
                              self._count_selection(dobj, [subset]),
                              cache = cache)

class AthenaPPGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        AMRGridPatch.__init__(self, id, filename = index.index_filename,
                              index = index)
        self.Parent = None
        self.Children = []
        self.Level = level

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
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
    _dataset_type='athena_pp'
    _data_file = None

    def __init__(self, ds, dataset_type='athena_pp'):
        self.dataset = weakref.proxy(ds)
        self.directory = os.path.dirname(self.dataset.filename)
        self.dataset_type = dataset_type
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.filename
        self._handle = ds._handle
        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        self.field_list = [("athena_pp", k) for k in self.dataset._field_map]

    def _count_grids(self):
        self.num_grids = self._handle.attrs["NumMeshBlocks"]

    def _parse_index(self):
        num_grids = self._handle.attrs["NumMeshBlocks"]

        self.grid_left_edge = np.zeros((num_grids, 3), dtype='float64')
        self.grid_right_edge = np.zeros((num_grids, 3), dtype='float64')
        self.grid_dimensions = np.zeros((num_grids, 3), dtype='int32')

        for i in range(num_grids):
            x = self._handle["x1f"][i,:]
            y = self._handle["x2f"][i,:]
            z = self._handle["x3f"][i,:]
            self.grid_left_edge[i] = np.array([x[0], y[0], z[0]], dtype='float64')
            self.grid_right_edge[i] = np.array([x[-1], y[-1], z[-1]], dtype='float64')
            self.grid_dimensions[i] = self._handle.attrs["MeshBlockSize"]
        levels = self._handle["Levels"][:]

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
        self.max_level = self._handle.attrs["MaxLevel"]

    def _chunk_io(self, dobj, cache = True, local_only = False):
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in gobjs:
            yield YTDataChunk(dobj, "io", [subset],
                              self._count_selection(dobj, [subset]),
                              cache = cache)

class AthenaPPDataset(Dataset):
    _field_info_class = AthenaPPFieldInfo
    _dataset_type = "athena_pp"

    def __init__(self, filename, dataset_type='athena_pp',
                 storage_filename=None, parameters=None,
                 units_override=None, unit_system="code"):
        self.fluid_types += ("athena_pp",)
        if parameters is None:
            parameters = {}
        self.specified_parameters = parameters
        if units_override is None:
            units_override = {}
        self._handle = HDF5FileHandler(filename)
        xrat = self._handle.attrs["RootGridX1"][2]
        yrat = self._handle.attrs["RootGridX2"][2]
        zrat = self._handle.attrs["RootGridX3"][2]
        if xrat != 1.0 or yrat != 1.0 or zrat != 1.0:
            self._index_class = AthenaPPLogarithmicIndex
            self.logarithmic = True
        else:
            self._index_class = AthenaPPHierarchy
            self.logarithmic = False
        Dataset.__init__(self, filename, dataset_type, units_override=units_override,
                         unit_system=unit_system)
        self.filename = filename
        if storage_filename is None:
            storage_filename = '%s.yt' % filename.split('/')[-1]
        self.storage_filename = storage_filename
        self.backup_filename = self.filename[:-4] + "_backup.gdf"

    def _set_code_unit_attributes(self):
        """
        Generates the conversion to various physical _units based on the
        parameter file
        """
        if "length_unit" not in self.units_override:
            self.no_cgs_equiv_length = True
        for unit, cgs in [("length", "cm"), ("time", "s"), ("mass", "g"),
                          ("temperature", "K")]:
            # We set these to cgs for now, but they may have been overridden
            if getattr(self, unit+'_unit', None) is not None:
                continue
            mylog.warning("Assuming 1.0 = 1.0 %s", cgs)
            setattr(self, "%s_unit" % unit, self.quan(1.0, cgs))

        self.magnetic_unit = np.sqrt(4*np.pi * self.mass_unit /
                                     (self.time_unit**2 * self.length_unit))
        self.magnetic_unit.convert_to_units("gauss")
        self.velocity_unit = self.length_unit / self.time_unit

    def _parse_parameter_file(self):

        xmin, xmax = self._handle.attrs["RootGridX1"][:2]
        ymin, ymax = self._handle.attrs["RootGridX2"][:2]
        zmin, zmax = self._handle.attrs["RootGridX3"][:2]

        self.domain_left_edge = np.array([xmin, ymin, zmin], dtype='float64')
        self.domain_right_edge = np.array([xmax, ymax, zmax], dtype='float64')

        self.geometry = geom_map[self._handle.attrs["Coordinates"].decode('utf-8')]
        self.domain_width = self.domain_right_edge-self.domain_left_edge
        self.domain_dimensions = self._handle.attrs["RootGridSize"]

        self._field_map = {}
        k = 0
        for (i, dname), num_var in zip(enumerate(self._handle.attrs["DatasetNames"]),
                                       self._handle.attrs["NumVariables"]):
            for j in range(num_var):
                fname = self._handle.attrs["VariableNames"][k].decode("ascii","ignore")
                self._field_map[fname] = (dname.decode("ascii","ignore"), j)
                k += 1

        self.refine_by = 2
        dimensionality = 3
        if self.domain_dimensions[2] == 1:
            dimensionality = 2
        if self.domain_dimensions[1] == 1:
            dimensionality = 1
        self.dimensionality = dimensionality
        self.current_time = self._handle.attrs["Time"]
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

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            if args[0].endswith('athdf'):
                return True
        except:
            pass
        return False

    @property
    def _skip_cache(self):
        return True

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]
