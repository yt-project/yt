import weakref
from typing import List, Tuple

import numpy as np

import yt.geometry.particle_deposit as particle_deposit
from yt._maintenance.deprecation import issue_deprecation_warning
from yt.config import ytcfg
from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer,
)
from yt.funcs import is_sequence
from yt.geometry.selection_routines import convert_mask_to_indices
from yt.units.yt_array import YTArray
from yt.utilities.exceptions import (
    YTFieldTypeNotFound,
    YTParticleDepositionNotImplemented,
)
from yt.utilities.lib.interpolators import ghost_zone_interpolate
from yt.utilities.lib.mesh_utilities import clamp_edges
from yt.utilities.nodal_data_utils import get_nodal_slices


class AMRGridPatch(YTSelectionContainer):
    _spatial = True
    _num_ghost_zones = 0
    _grids = None
    _id_offset = 1
    _cache_mask = True

    _type_name = "grid"
    _skip_add = True
    _con_args = ("id", "filename")
    _container_fields = (
        ("index", "dx"),
        ("index", "dy"),
        ("index", "dz"),
        ("index", "x"),
        ("index", "y"),
        ("index", "z"),
    )
    OverlappingSiblings = None

    def __init__(self, id, filename=None, index=None):
        super().__init__(index.dataset, None)
        self.id = id
        self._child_mask = self._child_indices = self._child_index_mask = None
        self.ds = index.dataset
        self._index = weakref.proxy(index)
        self.start_index = None
        self.filename = filename
        self._last_mask = None
        self._last_count = -1
        self._last_selector_id = None

    def get_global_startindex(self):
        """
        Return the integer starting index for each dimension at the current
        level.

        """
        if self.start_index is not None:
            return self.start_index
        if self.Parent is None:
            left = self.LeftEdge.d - self.ds.domain_left_edge.d
            start_index = left / self.dds.d
            return np.rint(start_index).astype("int64").ravel()

        pdx = self.Parent.dds.d
        di = np.rint((self.LeftEdge.d - self.Parent.LeftEdge.d) / pdx)
        start_index = self.Parent.get_global_startindex() + di
        self.start_index = (start_index * self.ds.refine_by).astype("int64").ravel()
        return self.start_index

    def __getitem__(self, key):
        tr = super().__getitem__(key)
        try:
            fields = self._determine_fields(key)
        except YTFieldTypeNotFound:
            return tr
        finfo = self.ds._get_field_info(*fields[0])
        if not finfo.sampling_type == "particle":
            num_nodes = 2 ** sum(finfo.nodal_flag)
            new_shape = list(self.ActiveDimensions)
            if num_nodes > 1:
                new_shape += [num_nodes]
            return tr.reshape(new_shape)
        return tr

    def convert(self, datatype):
        """
        This will attempt to convert a given unit to cgs from code units. It
        either returns the multiplicative factor or throws a KeyError.

        """
        return self.ds[datatype]

    @property
    def shape(self):
        return self.ActiveDimensions

    def _reshape_vals(self, arr):
        if len(arr.shape) == 3:
            return arr
        return arr.reshape(self.ActiveDimensions, order="C")

    def _generate_container_field(self, field):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        if field == ("index", "dx"):
            tr = self._current_chunk.fwidth[:, 0]
        elif field == ("index", "dy"):
            tr = self._current_chunk.fwidth[:, 1]
        elif field == ("index", "dz"):
            tr = self._current_chunk.fwidth[:, 2]
        elif field == ("index", "x"):
            tr = self._current_chunk.fcoords[:, 0]
        elif field == ("index", "y"):
            tr = self._current_chunk.fcoords[:, 1]
        elif field == ("index", "z"):
            tr = self._current_chunk.fcoords[:, 2]
        return self._reshape_vals(tr)

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz, at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        ds = self.ds
        index = self.index
        if self.Parent is not None:
            if not hasattr(self.Parent, "dds"):
                self.Parent._setup_dx()
            self.dds = self.Parent.dds.d / self.ds.refine_by
        else:
            LE, RE = (index.grid_left_edge[id, :].d, index.grid_right_edge[id, :].d)
            self.dds = (RE - LE) / self.ActiveDimensions
        if self.ds.dimensionality < 3:
            self.dds[2] = ds.domain_right_edge[2] - ds.domain_left_edge[2]
        elif self.ds.dimensionality < 2:
            self.dds[1] = ds.domain_right_edge[1] - ds.domain_left_edge[1]
        self.dds = self.dds.view(YTArray)
        self.dds.units = self.index.grid_left_edge.units

    def __repr__(self):
        cls_name = self.__class__.__name__
        return f"{cls_name}_{self.id:04d} ({self.ActiveDimensions})"

    def __int__(self):
        return self.id

    def clear_data(self):
        """
        Clear out the following things: child_mask, child_indices, all fields,
        all field parameters.

        """
        super().clear_data()
        self._setup_dx()

    def _prepare_grid(self):
        """Copies all the appropriate attributes from the index."""
        # This is definitely the slowest part of generating the index
        # Now we give it pointers to all of its attributes
        # Note that to keep in line with Enzo, we have broken PEP-8
        h = self.index  # cache it
        my_ind = self.id - self._id_offset
        self.ActiveDimensions = h.grid_dimensions[my_ind]
        self.LeftEdge = h.grid_left_edge[my_ind]
        self.RightEdge = h.grid_right_edge[my_ind]
        # This can be expensive so we allow people to disable this behavior
        # via a config option
        if ytcfg.get("yt", "reconstruct_index"):
            if is_sequence(self.Parent) and len(self.Parent) > 0:
                p = self.Parent[0]
            else:
                p = self.Parent
            if p is not None and p != []:
                # clamp grid edges to an integer multiple of the parent cell
                # width
                clamp_edges(self.LeftEdge, p.LeftEdge, p.dds)
                clamp_edges(self.RightEdge, p.RightEdge, p.dds)
        h.grid_levels[my_ind, 0] = self.Level
        # This might be needed for streaming formats
        # self.Time = h.gridTimes[my_ind,0]
        self.NumberOfParticles = h.grid_particle_count[my_ind, 0]

    def get_position(self, index):
        """Returns center position of an *index*."""
        pos = (index + 0.5) * self.dds + self.LeftEdge
        return pos

    def _fill_child_mask(self, child, mask, tofill, dlevel=1):
        rf = self.ds.refine_by
        if dlevel != 1:
            rf = rf**dlevel
        gi, cgi = self.get_global_startindex(), child.get_global_startindex()
        startIndex = np.maximum(0, cgi // rf - gi)
        endIndex = np.minimum(
            (cgi + child.ActiveDimensions) // rf - gi, self.ActiveDimensions
        )
        endIndex += startIndex == endIndex
        mask[
            startIndex[0] : endIndex[0],
            startIndex[1] : endIndex[1],
            startIndex[2] : endIndex[2],
        ] = tofill

    @property
    def child_mask(self):
        """
        Generates self.child_mask, which is zero where child grids exist (and
        thus, where higher resolution data is available).

        """
        child_mask = np.ones(self.ActiveDimensions, "bool")
        for child in self.Children:
            self._fill_child_mask(child, child_mask, 0)
        for sibling in self.OverlappingSiblings or []:
            self._fill_child_mask(sibling, child_mask, 0, dlevel=0)
        return child_mask

    @property
    def child_indices(self):
        return self.child_mask == 0

    @property
    def child_index_mask(self):
        """
        Generates self.child_index_mask, which is -1 where there is no child,
        and otherwise has the ID of the grid that resides there.

        """
        child_index_mask = np.zeros(self.ActiveDimensions, "int32") - 1
        for child in self.Children:
            self._fill_child_mask(child, child_index_mask, child.id)
        for sibling in self.OverlappingSiblings or []:
            self._fill_child_mask(sibling, child_index_mask, sibling.id, dlevel=0)
        return child_index_mask

    def retrieve_ghost_zones(self, n_zones, fields, all_levels=False, smoothed=False):
        # We will attempt this by creating a datacube that is exactly bigger
        # than the grid by nZones*dx in each direction
        nl = self.get_global_startindex() - n_zones
        new_left_edge = nl * self.dds + self.ds.domain_left_edge

        # Something different needs to be done for the root grid, though
        level = self.Level
        if all_levels:
            level = self.index.max_level + 1
        kwargs = {
            "dims": self.ActiveDimensions + 2 * n_zones,
            "num_ghost_zones": n_zones,
            "use_pbar": False,
            "fields": fields,
        }
        # This should update the arguments to set the field parameters to be
        # those of this grid.
        field_parameters = {}
        field_parameters.update(self.field_parameters)
        if smoothed:
            cube = self.ds.smoothed_covering_grid(
                level, new_left_edge, field_parameters=field_parameters, **kwargs
            )
        else:
            cube = self.ds.covering_grid(
                level, new_left_edge, field_parameters=field_parameters, **kwargs
            )
        cube._base_grid = self
        return cube

    def get_vertex_centered_data(
        self,
        fields: List[Tuple[str, str]],
        smoothed: bool = True,
        no_ghost: bool = False,
    ):
        _old_api = isinstance(fields, (str, tuple))
        if _old_api:
            issue_deprecation_warning(
                "get_vertex_centered_data() requires list of fields, rather than "
                "a single field as an argument.",
                since="3.3",
                removal="4.2",
            )
            fields = [fields]  # type: ignore

        # Make sure the field list has only unique entries
        fields = list(set(fields))
        new_fields = {}
        for field in fields:
            finfo = self.ds._get_field_info(field)
            new_fields[field] = self.ds.arr(
                np.zeros(self.ActiveDimensions + 1), finfo.output_units
            )
        if no_ghost:
            for field in fields:
                # Ensure we have the native endianness in this array.  Avoid making
                # a copy if possible.
                old_field = np.asarray(self[field], dtype="=f8")
                # We'll use the ghost zone routine, which will naturally
                # extrapolate here.
                input_left = np.array([0.5, 0.5, 0.5], dtype="float64")
                output_left = np.array([0.0, 0.0, 0.0], dtype="float64")
                # rf = 1 here
                ghost_zone_interpolate(
                    1, old_field, input_left, new_fields[field], output_left
                )
        else:
            cg = self.retrieve_ghost_zones(1, fields, smoothed=smoothed)
            for field in fields:
                src = cg[field].in_units(new_fields[field].units).d
                dest = new_fields[field].d
                np.add(dest, src[1:, 1:, 1:], dest)
                np.add(dest, src[:-1, 1:, 1:], dest)
                np.add(dest, src[1:, :-1, 1:], dest)
                np.add(dest, src[1:, 1:, :-1], dest)
                np.add(dest, src[:-1, 1:, :-1], dest)
                np.add(dest, src[1:, :-1, :-1], dest)
                np.add(dest, src[:-1, :-1, 1:], dest)
                np.add(dest, src[:-1, :-1, :-1], dest)
                np.multiply(dest, 0.125, dest)

        if _old_api:
            return new_fields[fields[0]]
        return new_fields

    def select_icoords(self, dobj):
        mask = self._get_selector_mask(dobj.selector)
        if mask is None:
            return np.empty((0, 3), dtype="int64")
        coords = convert_mask_to_indices(mask, self._last_count)
        coords += self.get_global_startindex()[None, :]
        return coords

    def select_fcoords(self, dobj):
        mask = self._get_selector_mask(dobj.selector)
        if mask is None:
            return np.empty((0, 3), dtype="float64")
        coords = convert_mask_to_indices(mask, self._last_count).astype("float64")
        coords += 0.5
        coords *= self.dds.d[None, :]
        coords += self.LeftEdge.d[None, :]
        return coords

    def select_fwidth(self, dobj):
        count = self.count(dobj.selector)
        if count == 0:
            return np.empty((0, 3), dtype="float64")
        coords = np.empty((count, 3), dtype="float64")
        for axis in range(3):
            coords[:, axis] = self.dds.d[axis]
        return coords

    def select_ires(self, dobj):
        mask = self._get_selector_mask(dobj.selector)
        if mask is None:
            return np.empty(0, dtype="int64")
        coords = np.empty(self._last_count, dtype="int64")
        coords[:] = self.Level
        return coords

    def select_tcoords(self, dobj):
        dt, t = dobj.selector.get_dt(self)
        return dt, t

    def smooth(self, *args, **kwargs):
        raise NotImplementedError

    def particle_operation(self, *args, **kwargs):
        raise NotImplementedError

    def deposit(self, positions, fields=None, method=None, kernel_name="cubic"):
        # Here we perform our particle deposition.
        cls = getattr(particle_deposit, f"deposit_{method}", None)
        if cls is None:
            raise YTParticleDepositionNotImplemented(method)
        # We allocate number of zones, not number of octs. Everything
        # inside this is Fortran ordered because of the ordering in the
        # octree deposit routines, so we reverse it here to match the
        # convention there
        nvals = tuple(self.ActiveDimensions[::-1])
        # append a dummy dimension because we are only depositing onto
        # one grid
        op = cls(nvals + (1,), kernel_name)
        op.initialize()
        if positions.size > 0:
            op.process_grid(self, positions, fields)
        vals = op.finalize()
        if vals is None:
            return
        # Fortran-ordered, so transpose.
        vals = vals.transpose()
        # squeeze dummy dimension we appended above
        return np.squeeze(vals, axis=0)

    def select_blocks(self, selector):
        mask = self._get_selector_mask(selector)
        yield self, mask

    def _get_selector_mask(self, selector):
        if self._cache_mask and hash(selector) == self._last_selector_id:
            mask = self._last_mask
        else:
            mask = selector.fill_mask_regular_grid(self)
            if self._cache_mask:
                self._last_mask = mask
            self._last_selector_id = hash(selector)
            if mask is None:
                self._last_count = 0
            else:
                self._last_count = mask.sum()
        return mask

    def select(self, selector, source, dest, offset):
        mask = self._get_selector_mask(selector)
        count = self.count(selector)
        if count == 0:
            return 0
        dim = np.squeeze(self.ds.dimensionality)
        nodal_flag = source.shape[:dim] - self.ActiveDimensions[:dim]
        if sum(nodal_flag) == 0:
            dest[offset : offset + count] = source[mask]
        else:
            slices = get_nodal_slices(source.shape, nodal_flag, dim)
            for i, sl in enumerate(slices):
                dest[offset : offset + count, i] = source[tuple(sl)][np.squeeze(mask)]
        return count

    def count(self, selector):
        mask = self._get_selector_mask(selector)
        if mask is None:
            return 0
        return self._last_count

    def count_particles(self, selector, x, y, z):
        # We don't cache the selector results
        count = selector.count_points(x, y, z, 0.0)
        return count

    def select_particles(self, selector, x, y, z):
        mask = selector.select_points(x, y, z, 0.0)
        return mask
