import os
import sys
import time
import uuid
import weakref
from collections import UserDict
from itertools import chain, product, repeat
from numbers import Number as numeric_type
from typing import Type

import numpy as np
from more_itertools import always_iterable

from yt.data_objects.field_data import YTFieldData
from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.index_subobjects.octree_subset import OctreeSubset
from yt.data_objects.index_subobjects.stretched_grid import StretchedGrid
from yt.data_objects.index_subobjects.unstructured_mesh import (
    SemiStructuredMesh,
    UnstructuredMesh,
)
from yt.data_objects.particle_unions import ParticleUnion
from yt.data_objects.static_output import Dataset, ParticleFile
from yt.data_objects.unions import MeshUnion
from yt.frontends.sph.data_structures import SPHParticleIndex
from yt.geometry.geometry_handler import Index, YTDataChunk
from yt.geometry.grid_geometry_handler import GridIndex
from yt.geometry.oct_container import OctreeContainer
from yt.geometry.oct_geometry_handler import OctreeIndex
from yt.geometry.unstructured_mesh_handler import UnstructuredIndex
from yt.units import YTQuantity
from yt.utilities.io_handler import io_registry
from yt.utilities.lib.cykdtree import PyKDTree
from yt.utilities.lib.misc_utilities import (
    _obtain_coords_and_widths,
    get_box_grids_level,
)
from yt.utilities.lib.particle_kdtree_tools import (
    estimate_density,
    generate_smoothing_length,
)
from yt.utilities.logger import ytLogger as mylog

from .definitions import process_data, set_particle_types
from .fields import StreamFieldInfo

if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    from yt._maintenance.backports import cached_property


class StreamGrid(AMRGridPatch):
    """
    Class representing a single In-memory Grid instance.
    """

    __slots__ = ["proc_num"]
    _id_offset = 0

    def __init__(self, id, index):
        """
        Returns an instance of StreamGrid with *id*, associated with *filename*
        and *index*.
        """
        # All of the field parameters will be passed to us as needed.
        AMRGridPatch.__init__(self, id, filename=None, index=index)
        self._children_ids = []
        self._parent_id = -1
        self.Level = -1

    def set_filename(self, filename):
        pass

    @property
    def Parent(self):
        if self._parent_id == -1:
            return None
        return self.index.grids[self._parent_id - self._id_offset]

    @property
    def Children(self):
        return [self.index.grids[cid - self._id_offset] for cid in self._children_ids]


class StreamStretchedGrid(StretchedGrid):
    _id_offset = 0

    def __init__(self, id, index):
        cell_widths = index.grid_cell_widths[id - self._id_offset]
        super().__init__(id, cell_widths, index=index)
        self._children_ids = []
        self._parent_id = -1
        self.Level = -1

    @property
    def Parent(self):
        if self._parent_id == -1:
            return None
        return self.index.grids[self._parent_id - self._id_offset]

    @property
    def Children(self):
        return [self.index.grids[cid - self._id_offset] for cid in self._children_ids]


class StreamHandler:
    def __init__(
        self,
        left_edges,
        right_edges,
        dimensions,
        levels,
        parent_ids,
        particle_count,
        processor_ids,
        fields,
        field_units,
        code_units,
        io=None,
        particle_types=None,
        periodicity=(True, True, True),
        *,
        cell_widths=None,
        parameters=None,
    ):
        if particle_types is None:
            particle_types = {}
        self.left_edges = np.array(left_edges)
        self.right_edges = np.array(right_edges)
        self.dimensions = dimensions
        self.levels = levels
        self.parent_ids = parent_ids
        self.particle_count = particle_count
        self.processor_ids = processor_ids
        self.num_grids = self.levels.size
        self.fields = fields
        self.field_units = field_units
        self.code_units = code_units
        self.io = io
        self.particle_types = particle_types
        self.periodicity = periodicity
        self.cell_widths = cell_widths

        if parameters is None:
            self.parameters = {}
        else:
            self.parameters = parameters.copy()

    def get_fields(self):
        return self.fields.all_fields

    def get_particle_type(self, field):

        if field in self.particle_types:
            return self.particle_types[field]
        else:
            return False


class StreamHierarchy(GridIndex):

    grid = StreamGrid

    def __init__(self, ds, dataset_type=None):
        self.dataset_type = dataset_type
        self.float_type = "float64"
        self.dataset = weakref.proxy(ds)  # for _obtain_enzo
        self.stream_handler = ds.stream_handler
        self.float_type = "float64"
        self.directory = os.getcwd()
        GridIndex.__init__(self, ds, dataset_type)

    def _count_grids(self):
        self.num_grids = self.stream_handler.num_grids

    def _icoords_to_fcoords(self, icoords, ires, axes=None):
        """
        We check here that we have cell_widths, and if we do, we will provide them.
        """
        if self.grid_cell_widths is None:
            return super()._icoords_to_fcoords(icoords, ires, axes)
        if axes is None:
            axes = [0, 1, 2]
        # Transpose these by reversing the shape
        coords = np.empty(icoords.shape, dtype="f8")
        cell_widths = np.empty(icoords.shape, dtype="f8")
        for i, ax in enumerate(axes):
            coords[:, i], cell_widths[:, i] = _obtain_coords_and_widths(
                icoords[:, i],
                ires,
                self.grid_cell_widths[0][ax],
                self.ds.domain_left_edge[ax].d,
            )
        return coords, cell_widths

    def _parse_index(self):
        self.grid_dimensions = self.stream_handler.dimensions
        self.grid_left_edge[:] = self.stream_handler.left_edges
        self.grid_right_edge[:] = self.stream_handler.right_edges
        self.grid_levels[:] = self.stream_handler.levels
        self.min_level = self.grid_levels.min()
        self.grid_procs = self.stream_handler.processor_ids
        self.grid_particle_count[:] = self.stream_handler.particle_count
        if self.stream_handler.cell_widths is not None:
            self.grid_cell_widths = self.stream_handler.cell_widths[:]
            self.grid = StreamStretchedGrid
        else:
            self.grid_cell_widths = None
        mylog.debug("Copying reverse tree")
        self.grids = []
        # We enumerate, so it's 0-indexed id and 1-indexed pid
        for id in range(self.num_grids):
            self.grids.append(self.grid(id, self))
            self.grids[id].Level = self.grid_levels[id, 0]
        parent_ids = self.stream_handler.parent_ids
        if parent_ids is not None:
            reverse_tree = self.stream_handler.parent_ids.tolist()
            # Initial setup:
            for gid, pid in enumerate(reverse_tree):
                if pid >= 0:
                    self.grids[gid]._parent_id = pid
                    self.grids[pid]._children_ids.append(self.grids[gid].id)
        else:
            mylog.debug("Reconstructing parent-child relationships")
            self._reconstruct_parent_child()
        self.max_level = self.grid_levels.max()
        mylog.debug("Preparing grids")
        temp_grids = np.empty(self.num_grids, dtype="object")
        for i, grid in enumerate(self.grids):
            if (i % 1e4) == 0:
                mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)
            grid.filename = None
            grid._prepare_grid()
            grid._setup_dx()
            grid.proc_num = self.grid_procs[i]
            temp_grids[i] = grid
        self.grids = temp_grids
        mylog.debug("Prepared")

    def _reconstruct_parent_child(self):
        mask = np.empty(len(self.grids), dtype="int32")
        mylog.debug("First pass; identifying child grids")
        for i, grid in enumerate(self.grids):
            get_box_grids_level(
                self.grid_left_edge[i, :],
                self.grid_right_edge[i, :],
                self.grid_levels[i] + 1,
                self.grid_left_edge,
                self.grid_right_edge,
                self.grid_levels,
                mask,
            )
            ids = np.where(mask.astype("bool"))
            grid._children_ids = ids[0]  # where is a tuple
        mylog.debug("Second pass; identifying parents")
        self.stream_handler.parent_ids = (
            np.zeros(self.stream_handler.num_grids, "int64") - 1
        )
        for i, grid in enumerate(self.grids):  # Second pass
            for child in grid.Children:
                child._parent_id = i
                # _id_offset = 0
                self.stream_handler.parent_ids[child.id] = i

    def _initialize_grid_arrays(self):
        GridIndex._initialize_grid_arrays(self)
        self.grid_procs = np.zeros((self.num_grids, 1), "int32")

    def _detect_output_fields(self):
        # NOTE: Because particle unions add to the actual field list, without
        # having the keys in the field list itself, we need to double check
        # here.
        fl = set(self.stream_handler.get_fields())
        fl.update(set(getattr(self, "field_list", [])))
        self.field_list = list(fl)

    def _populate_grid_objects(self):
        for g in self.grids:
            g._setup_dx()
        self.max_level = self.grid_levels.max()

    def _setup_data_io(self):
        if self.stream_handler.io is not None:
            self.io = self.stream_handler.io
        else:
            self.io = io_registry[self.dataset_type](self.ds)

    def _reset_particle_count(self):
        self.grid_particle_count[:] = self.stream_handler.particle_count
        for i, grid in enumerate(self.grids):
            grid.NumberOfParticles = self.grid_particle_count[i, 0]

    def update_data(self, data):
        """
        Update the stream data with a new data dict. If fields already exist,
        they will be replaced, but if they do not, they will be added. Fields
        already in the stream but not part of the data dict will be left
        alone.
        """
        particle_types = set_particle_types(data[0])

        self.stream_handler.particle_types.update(particle_types)
        self.ds._find_particle_types()

        for i, grid in enumerate(self.grids):
            field_units, gdata, number_of_particles = process_data(data[i])
            self.stream_handler.particle_count[i] = number_of_particles
            self.stream_handler.field_units.update(field_units)
            for field in gdata:
                if field in grid.field_data:
                    grid.field_data.pop(field, None)
                self.stream_handler.fields[grid.id][field] = gdata[field]

        self._reset_particle_count()
        # We only want to create a superset of fields here.
        for field in self.ds.field_list:
            if field[0] == "all":
                self.ds.field_list.remove(field)
        self._detect_output_fields()
        self.ds.create_field_info()
        mylog.debug("Creating Particle Union 'all'")
        pu = ParticleUnion("all", list(self.ds.particle_types_raw))
        self.ds.add_particle_union(pu)
        self.ds.particle_types = tuple(set(self.ds.particle_types))


class StreamDataset(Dataset):
    _index_class: Type[Index] = StreamHierarchy
    _field_info_class = StreamFieldInfo
    _dataset_type = "stream"

    def __init__(
        self,
        stream_handler,
        storage_filename=None,
        geometry="cartesian",
        unit_system="cgs",
        default_species_fields=None,
    ):
        self.fluid_types += ("stream",)
        self.geometry = geometry
        self.stream_handler = stream_handler
        self._find_particle_types()
        name = f"InMemoryParameterFile_{uuid.uuid4().hex}"
        from yt.data_objects.static_output import _cached_datasets

        _cached_datasets[name] = self
        Dataset.__init__(
            self,
            name,
            self._dataset_type,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
        )

    @property
    def filename(self):
        return self.stream_handler.name

    @cached_property
    def unique_identifier(self) -> str:
        return str(self.parameters["CurrentTimeIdentifier"])

    def _parse_parameter_file(self):
        self.parameters["CurrentTimeIdentifier"] = time.time()
        self.domain_left_edge = self.stream_handler.domain_left_edge.copy()
        self.domain_right_edge = self.stream_handler.domain_right_edge.copy()
        self.refine_by = self.stream_handler.refine_by
        self.dimensionality = self.stream_handler.dimensionality
        self._periodicity = self.stream_handler.periodicity
        self.domain_dimensions = self.stream_handler.domain_dimensions
        self.current_time = self.stream_handler.simulation_time
        self.gamma = 5.0 / 3.0
        self.parameters["EOSType"] = -1
        self.parameters["CosmologyHubbleConstantNow"] = 1.0
        self.parameters["CosmologyCurrentRedshift"] = 1.0
        self.parameters["HydroMethod"] = -1
        self.parameters.update(self.stream_handler.parameters)
        if self.stream_handler.cosmology_simulation:
            self.cosmological_simulation = 1
            self.current_redshift = self.stream_handler.current_redshift
            self.omega_lambda = self.stream_handler.omega_lambda
            self.omega_matter = self.stream_handler.omega_matter
            self.hubble_constant = self.stream_handler.hubble_constant
        else:
            self.current_redshift = 0.0
            self.omega_lambda = 0.0
            self.omega_matter = 0.0
            self.hubble_constant = 0.0
            self.cosmological_simulation = 0

    def _set_units(self):
        self.field_units = self.stream_handler.field_units

    def _set_code_unit_attributes(self):
        base_units = self.stream_handler.code_units
        attrs = (
            "length_unit",
            "mass_unit",
            "time_unit",
            "velocity_unit",
            "magnetic_unit",
        )
        cgs_units = ("cm", "g", "s", "cm/s", "gauss")
        for unit, attr, cgs_unit in zip(base_units, attrs, cgs_units):
            if isinstance(unit, str):
                if unit == "code_magnetic":
                    # If no magnetic unit was explicitly specified
                    # we skip it now and take care of it at the bottom
                    continue
                else:
                    uq = self.quan(1.0, unit)
            elif isinstance(unit, numeric_type):
                uq = self.quan(unit, cgs_unit)
            elif isinstance(unit, YTQuantity):
                uq = unit
            elif isinstance(unit, tuple):
                uq = self.quan(unit[0], unit[1])
            else:
                raise RuntimeError(f"{attr} ({unit}) is invalid.")
            setattr(self, attr, uq)
        if not hasattr(self, "magnetic_unit"):
            self.magnetic_unit = np.sqrt(
                4 * np.pi * self.mass_unit / (self.time_unit**2 * self.length_unit)
            )

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        return False

    @property
    def _skip_cache(self):
        return True

    def _find_particle_types(self):
        particle_types = set()
        for k, v in self.stream_handler.particle_types.items():
            if v:
                particle_types.add(k[0])
        self.particle_types = tuple(particle_types)
        self.particle_types_raw = self.particle_types


class StreamDictFieldHandler(UserDict):
    _additional_fields = ()

    @property
    def all_fields(self):
        self_fields = chain.from_iterable(s.keys() for s in self.values())
        self_fields = list(set(self_fields))
        fields = list(self._additional_fields) + self_fields
        fields = list(set(fields))
        return fields


class StreamParticleIndex(SPHParticleIndex):
    def __init__(self, ds, dataset_type=None):
        self.stream_handler = ds.stream_handler
        super().__init__(ds, dataset_type)

    def _setup_data_io(self):
        if self.stream_handler.io is not None:
            self.io = self.stream_handler.io
        else:
            self.io = io_registry[self.dataset_type](self.ds)

    def update_data(self, data):
        """
        Update the stream data with a new data dict. If fields already exist,
        they will be replaced, but if they do not, they will be added. Fields
        already in the stream but not part of the data dict will be left
        alone.
        """
        # Alias
        ds = self.ds
        handler = ds.stream_handler

        # Preprocess
        field_units, data, _ = process_data(data)
        pdata = {}
        for key in data.keys():
            if not isinstance(key, tuple):
                field = ("io", key)
                mylog.debug("Reassigning '%s' to '%s'", key, field)
            else:
                field = key
            pdata[field] = data[key]
        data = pdata  # Drop reference count
        particle_types = set_particle_types(data)

        # Update particle types
        handler.particle_types.update(particle_types)
        ds._find_particle_types()

        # Update fields
        handler.field_units.update(field_units)
        fields = handler.fields
        for field in data.keys():
            if field not in fields._additional_fields:
                fields._additional_fields += (field,)
        fields["stream_file"].update(data)

        # Update field list
        for field in self.ds.field_list:
            if field[0] in ["all", "nbody"]:
                self.ds.field_list.remove(field)
        self._detect_output_fields()
        self.ds.create_field_info()


class StreamParticleFile(ParticleFile):
    pass


class StreamParticlesDataset(StreamDataset):
    _index_class = StreamParticleIndex
    _file_class = StreamParticleFile
    _field_info_class = StreamFieldInfo
    _dataset_type = "stream_particles"
    file_count = 1
    filename_template = "stream_file"
    _proj_type = "particle_proj"

    def __init__(
        self,
        stream_handler,
        storage_filename=None,
        geometry="cartesian",
        unit_system="cgs",
        default_species_fields=None,
    ):
        super().__init__(
            stream_handler,
            storage_filename=storage_filename,
            geometry=geometry,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
        )
        fields = list(stream_handler.fields["stream_file"].keys())
        # This is the current method of detecting SPH data.
        # This should be made more flexible in the future.
        if ("io", "density") in fields and ("io", "smoothing_length") in fields:
            self._sph_ptypes = ("io",)

    def add_sph_fields(self, n_neighbors=32, kernel="cubic", sph_ptype="io"):
        """Add SPH fields for the specified particle type.

        For a particle type with "particle_position" and "particle_mass" already
        defined, this method adds the "smoothing_length" and "density" fields.
        "smoothing_length" is computed as the distance to the nth nearest
        neighbor. "density" is computed as the SPH (gather) smoothed mass. The
        SPH fields are added only if they don't already exist.

        Parameters
        ----------
        n_neighbors : int
            The number of neighbors to use in smoothing length computation.
        kernel : str
            The kernel function to use in density estimation.
        sph_ptype : str
            The SPH particle type. Each dataset has one sph_ptype only. This
            method will overwrite existing sph_ptype of the dataset.

        """
        mylog.info("Generating SPH fields")

        # Unify units
        l_unit = "code_length"
        m_unit = "code_mass"
        d_unit = "code_mass / code_length**3"

        # Read basic fields
        ad = self.all_data()
        pos = ad[sph_ptype, "particle_position"].to(l_unit).d
        mass = ad[sph_ptype, "particle_mass"].to(m_unit).d

        # Construct k-d tree
        kdtree = PyKDTree(
            pos.astype("float64"),
            left_edge=self.domain_left_edge.to_value(l_unit),
            right_edge=self.domain_right_edge.to_value(l_unit),
            periodic=self.periodicity,
            leafsize=2 * int(n_neighbors),
        )
        order = np.argsort(kdtree.idx)

        def exists(fname):
            if (sph_ptype, fname) in self.derived_field_list:
                mylog.info(
                    "Field ('%s','%s') already exists. Skipping", sph_ptype, fname
                )
                return True
            else:
                mylog.info("Generating field ('%s','%s')", sph_ptype, fname)
                return False

        data = {}

        # Add smoothing length field
        fname = "smoothing_length"
        if not exists(fname):
            hsml = generate_smoothing_length(pos[kdtree.idx], kdtree, n_neighbors)
            hsml = hsml[order]
            data[(sph_ptype, "smoothing_length")] = (hsml, l_unit)
        else:
            hsml = ad[sph_ptype, fname].to(l_unit).d

        # Add density field
        fname = "density"
        if not exists(fname):
            dens = estimate_density(
                pos[kdtree.idx],
                mass[kdtree.idx],
                hsml[kdtree.idx],
                kdtree,
                kernel_name=kernel,
            )
            dens = dens[order]
            data[(sph_ptype, "density")] = (dens, d_unit)

        # Add fields
        self._sph_ptypes = (sph_ptype,)
        self.index.update_data(data)
        self.num_neighbors = n_neighbors


_cis = np.fromiter(
    chain.from_iterable(product([0, 1], [0, 1], [0, 1])), dtype=np.int64, count=8 * 3
)
_cis.shape = (8, 3)


def hexahedral_connectivity(xgrid, ygrid, zgrid):
    r"""Define the cell coordinates and cell neighbors of a hexahedral mesh
    for a semistructured grid. Used to specify the connectivity and
    coordinates parameters used in
    :func:`~yt.frontends.stream.data_structures.load_hexahedral_mesh`.

    Parameters
    ----------
    xgrid : array_like
       x-coordinates of boundaries of the hexahedral cells. Should be a
       one-dimensional array.
    ygrid : array_like
       y-coordinates of boundaries of the hexahedral cells. Should be a
       one-dimensional array.
    zgrid : array_like
       z-coordinates of boundaries of the hexahedral cells. Should be a
       one-dimensional array.

    Returns
    -------
    coords : array_like
        The list of (x,y,z) coordinates of the vertices of the mesh.
        Is of size (M,3) where M is the number of vertices.
    connectivity : array_like
        For each hexahedron h in the mesh, gives the index of each of h's
        neighbors. Is of size (N,8), where N is the number of hexahedra.

    Examples
    --------

    >>> xgrid = np.array([-1, -0.25, 0, 0.25, 1])
    >>> coords, conn = hexahedral_connectivity(xgrid, xgrid, xgrid)
    >>> coords
    array([[-1.  , -1.  , -1.  ],
           [-1.  , -1.  , -0.25],
           [-1.  , -1.  ,  0.  ],
           ...,
           [ 1.  ,  1.  ,  0.  ],
           [ 1.  ,  1.  ,  0.25],
           [ 1.  ,  1.  ,  1.  ]])

    >>> conn
    array([[  0,   1,   5,   6,  25,  26,  30,  31],
           [  1,   2,   6,   7,  26,  27,  31,  32],
           [  2,   3,   7,   8,  27,  28,  32,  33],
           ...,
           [ 91,  92,  96,  97, 116, 117, 121, 122],
           [ 92,  93,  97,  98, 117, 118, 122, 123],
           [ 93,  94,  98,  99, 118, 119, 123, 124]])
    """
    nx = len(xgrid)
    ny = len(ygrid)
    nz = len(zgrid)
    coords = np.zeros((nx, ny, nz, 3), dtype="float64", order="C")
    coords[:, :, :, 0] = xgrid[:, None, None]
    coords[:, :, :, 1] = ygrid[None, :, None]
    coords[:, :, :, 2] = zgrid[None, None, :]
    coords.shape = (nx * ny * nz, 3)
    cycle = np.rollaxis(np.indices((nx - 1, ny - 1, nz - 1)), 0, 4)
    cycle.shape = ((nx - 1) * (ny - 1) * (nz - 1), 3)
    off = _cis + cycle[:, np.newaxis]
    connectivity = np.array(
        ((off[:, :, 0] * ny) + off[:, :, 1]) * nz + off[:, :, 2], order="C"
    )
    return coords, connectivity


class StreamHexahedralMesh(SemiStructuredMesh):
    _connectivity_length = 8
    _index_offset = 0


class StreamHexahedralHierarchy(UnstructuredIndex):
    def __init__(self, ds, dataset_type=None):
        self.stream_handler = ds.stream_handler
        super().__init__(ds, dataset_type)

    def _initialize_mesh(self):
        coords = self.stream_handler.fields.pop("coordinates")
        connect = self.stream_handler.fields.pop("connectivity")
        self.meshes = [
            StreamHexahedralMesh(0, self.index_filename, connect, coords, self)
        ]

    def _setup_data_io(self):
        if self.stream_handler.io is not None:
            self.io = self.stream_handler.io
        else:
            self.io = io_registry[self.dataset_type](self.ds)

    def _detect_output_fields(self):
        self.field_list = list(set(self.stream_handler.get_fields()))


class StreamHexahedralDataset(StreamDataset):
    _index_class = StreamHexahedralHierarchy
    _field_info_class = StreamFieldInfo
    _dataset_type = "stream_hexahedral"


class StreamOctreeSubset(OctreeSubset):
    domain_id = 1
    _domain_offset = 1

    def __init__(self, base_region, ds, oct_handler, num_zones=2, num_ghost_zones=0):
        self._num_zones = num_zones
        self.field_data = YTFieldData()
        self.field_parameters = {}
        self.ds = ds
        self.oct_handler = oct_handler
        self._last_mask = None
        self._last_selector_id = None
        self._current_particle_type = "io"
        self._current_fluid_type = self.ds.default_fluid_type
        self.base_region = base_region
        self.base_selector = base_region.selector

        self._num_ghost_zones = num_ghost_zones

        if num_ghost_zones > 0:
            if not all(ds.periodicity):
                mylog.warning(
                    "Ghost zones will wrongly assume the domain to be periodic."
                )
            base_grid = StreamOctreeSubset(base_region, ds, oct_handler, num_zones)
            self._base_grid = base_grid

    def retrieve_ghost_zones(self, ngz, fields, smoothed=False):
        try:
            new_subset = self._subset_with_gz
            mylog.debug("Reusing previous subset with ghost zone.")
        except AttributeError:
            new_subset = StreamOctreeSubset(
                self.base_region,
                self.ds,
                self.oct_handler,
                self._num_zones,
                num_ghost_zones=ngz,
            )
            self._subset_with_gz = new_subset

        return new_subset

    def _fill_no_ghostzones(self, content, dest, selector, offset):
        # Here we get a copy of the file, which we skip through and read the
        # bits we want.
        oct_handler = self.oct_handler
        cell_count = selector.count_oct_cells(self.oct_handler, self.domain_id)
        levels, cell_inds, file_inds = self.oct_handler.file_index_octs(
            selector, self.domain_id, cell_count
        )
        levels[:] = 0
        dest.update((field, np.empty(cell_count, dtype="float64")) for field in content)
        # Make references ...
        count = oct_handler.fill_level(
            0, levels, cell_inds, file_inds, dest, content, offset
        )
        return count

    def _fill_with_ghostzones(self, content, dest, selector, offset):
        oct_handler = self.oct_handler
        ndim = self.ds.dimensionality
        cell_count = (
            selector.count_octs(self.oct_handler, self.domain_id) * self.nz**ndim
        )

        gz_cache = getattr(self, "_ghost_zone_cache", None)
        if gz_cache:
            levels, cell_inds, file_inds, domains = gz_cache
        else:
            gz_cache = (
                levels,
                cell_inds,
                file_inds,
                domains,
            ) = oct_handler.file_index_octs_with_ghost_zones(
                selector, self.domain_id, cell_count
            )
            self._ghost_zone_cache = gz_cache
        levels[:] = 0
        dest.update((field, np.empty(cell_count, dtype="float64")) for field in content)
        # Make references ...
        oct_handler.fill_level(0, levels, cell_inds, file_inds, dest, content, offset)

    def fill(self, content, dest, selector, offset):
        if self._num_ghost_zones == 0:
            return self._fill_no_ghostzones(content, dest, selector, offset)
        else:
            return self._fill_with_ghostzones(content, dest, selector, offset)


class StreamOctreeHandler(OctreeIndex):
    def __init__(self, ds, dataset_type=None):
        self.stream_handler = ds.stream_handler
        self.dataset_type = dataset_type
        super().__init__(ds, dataset_type)

    def _setup_data_io(self):
        if self.stream_handler.io is not None:
            self.io = self.stream_handler.io
        else:
            self.io = io_registry[self.dataset_type](self.ds)

    def _initialize_oct_handler(self):
        header = dict(
            dims=[1, 1, 1],
            left_edge=self.ds.domain_left_edge,
            right_edge=self.ds.domain_right_edge,
            octree=self.ds.octree_mask,
            num_zones=self.ds.num_zones,
            partial_coverage=self.ds.partial_coverage,
        )
        self.oct_handler = OctreeContainer.load_octree(header)

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_chunk_info", None) is None:
            base_region = getattr(dobj, "base_region", dobj)
            subset = [
                StreamOctreeSubset(
                    base_region,
                    self.dataset,
                    self.oct_handler,
                    self.ds.num_zones,
                )
            ]
            dobj._chunk_info = subset
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, None)

    def _chunk_spatial(self, dobj, ngz, sort=None, preload_fields=None):
        sobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        # This is where we will perform cutting of the Octree and
        # load-balancing.  That may require a specialized selector object to
        # cut based on some space-filling curve index.
        for og in sobjs:
            if ngz > 0:
                g = og.retrieve_ghost_zones(ngz, [], smoothed=True)
            else:
                g = og
            yield YTDataChunk(dobj, "spatial", [g])

    def _chunk_io(self, dobj, cache=True, local_only=False):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in oobjs:
            yield YTDataChunk(dobj, "io", [subset], None, cache=cache)

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        super()._setup_classes(dd)

    def _detect_output_fields(self):
        # NOTE: Because particle unions add to the actual field list, without
        # having the keys in the field list itself, we need to double check
        # here.
        fl = set(self.stream_handler.get_fields())
        fl.update(set(getattr(self, "field_list", [])))
        self.field_list = list(fl)


class StreamOctreeDataset(StreamDataset):
    _index_class = StreamOctreeHandler
    _field_info_class = StreamFieldInfo
    _dataset_type = "stream_octree"

    levelmax = None

    def __init__(
        self,
        stream_handler,
        storage_filename=None,
        geometry="cartesian",
        unit_system="cgs",
        default_species_fields=None,
    ):
        super().__init__(
            stream_handler,
            storage_filename,
            geometry,
            unit_system,
            default_species_fields=default_species_fields,
        )
        # Set up levelmax
        self.max_level = stream_handler.levels.max()
        self.min_level = stream_handler.levels.min()


class StreamUnstructuredMesh(UnstructuredMesh):
    _index_offset = 0

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._connectivity_length = self.connectivity_indices.shape[1]


class StreamUnstructuredIndex(UnstructuredIndex):
    def __init__(self, ds, dataset_type=None):
        self.stream_handler = ds.stream_handler
        super().__init__(ds, dataset_type)

    def _initialize_mesh(self):
        coords = self.stream_handler.fields.pop("coordinates")
        connect = always_iterable(self.stream_handler.fields.pop("connectivity"))

        self.meshes = [
            StreamUnstructuredMesh(i, self.index_filename, c1, c2, self)
            for i, (c1, c2) in enumerate(zip(connect, repeat(coords)))
        ]
        self.mesh_union = MeshUnion("mesh_union", self.meshes)

    def _setup_data_io(self):
        if self.stream_handler.io is not None:
            self.io = self.stream_handler.io
        else:
            self.io = io_registry[self.dataset_type](self.ds)

    def _detect_output_fields(self):
        self.field_list = list(set(self.stream_handler.get_fields()))
        fnames = list({fn for ft, fn in self.field_list})
        self.field_list += [("all", fname) for fname in fnames]


class StreamUnstructuredMeshDataset(StreamDataset):
    _index_class = StreamUnstructuredIndex
    _field_info_class = StreamFieldInfo
    _dataset_type = "stream_unstructured"

    def _find_particle_types(self):
        pass
