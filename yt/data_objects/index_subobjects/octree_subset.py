import sys
from contextlib import contextmanager
from itertools import product, repeat
from typing import Tuple

import numpy as np
from unyt import unyt_array

import yt.geometry.particle_deposit as particle_deposit
import yt.geometry.particle_smooth as particle_smooth
from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer,
)
from yt.geometry.particle_oct_container import ParticleOctreeContainer
from yt.units.dimensions import length  # type: ignore
from yt.utilities.exceptions import (
    YTFieldTypeNotFound,
    YTInvalidPositionArray,
    YTParticleDepositionNotImplemented,
)
from yt.utilities.lib.geometry_utils import compute_morton
from yt.utilities.logger import ytLogger as mylog

if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    from yt._maintenance.backports import cached_property


def cell_count_cache(func):
    def cc_cache_func(self, dobj):
        if hash(dobj.selector) != self._last_selector_id:
            self._cell_count = -1
        rv = func(self, dobj)
        self._cell_count = rv.shape[0]
        self._last_selector_id = hash(dobj.selector)
        return rv

    return cc_cache_func


class OctreeSubset(YTSelectionContainer):
    _spatial = True
    _num_ghost_zones = 0
    _type_name = "octree_subset"
    _skip_add = True
    _con_args: Tuple[str, ...] = ("base_region", "domain", "ds")
    _domain_offset = 0
    _cell_count = -1
    _block_order = "C"

    def __init__(self, base_region, domain, ds, num_zones=2, num_ghost_zones=0):
        super().__init__(ds, None)
        self._num_zones = num_zones
        self._num_ghost_zones = num_ghost_zones
        self.domain = domain
        self.domain_id = domain.domain_id
        self.ds = domain.ds
        self._index = self.ds.index
        self.oct_handler = domain.oct_handler
        self._last_mask = None
        self._last_selector_id = None
        self.base_region = base_region
        self.base_selector = base_region.selector

    def __getitem__(self, key):
        tr = super().__getitem__(key)
        try:
            fields = self._determine_fields(key)
        except YTFieldTypeNotFound:
            return tr
        finfo = self.ds._get_field_info(*fields[0])
        if not finfo.sampling_type == "particle":
            # We may need to reshape the field, if it is being queried from
            # field_data.  If it's already cached, it just passes through.
            if len(tr.shape) < 4:
                tr = self._reshape_vals(tr)
            return tr
        return tr

    @property
    def nz(self):
        return self._num_zones + 2 * self._num_ghost_zones

    def _reshape_vals(self, arr):
        nz = self.nz
        if len(arr.shape) <= 2:
            n_oct = arr.shape[0] // (nz**3)
        elif arr.shape[-1] == 3:
            n_oct = arr.shape[-2]
        else:
            n_oct = arr.shape[-1]
        if arr.size == nz * nz * nz * n_oct:
            new_shape = (nz, nz, nz, n_oct)
        elif arr.size == nz * nz * nz * n_oct * 3:
            new_shape = (nz, nz, nz, n_oct, 3)
        else:
            raise RuntimeError
        # Note that if arr is already F-contiguous, this *shouldn't* copy the
        # data.  But, it might.  However, I can't seem to figure out how to
        # make the assignment to .shape, which *won't* copy the data, make the
        # resultant array viewed in Fortran order.
        arr = arr.reshape(new_shape, order="F")
        return arr

    def mask_refinement(self, selector):
        mask = self.oct_handler.mask(selector, domain_id=self.domain_id)
        return mask

    def select_blocks(self, selector):
        mask = self.oct_handler.mask(selector, domain_id=self.domain_id)
        slicer = OctreeSubsetBlockSlice(self, self.ds)
        for i, sl in slicer:
            yield sl, np.atleast_3d(mask[i, ...])

    def select_tcoords(self, dobj):
        # These will not be pre-allocated, which can be a problem for speed and
        # memory usage.
        dts, ts = [], []
        for sl, mask in self.select_blocks(dobj.selector):
            sl.child_mask = np.asfortranarray(mask)
            dt, t = dobj.selector.get_dt(sl)
            dts.append(dt)
            ts.append(t)
        if len(dts) == len(ts) == 0:
            return np.empty(0, "f8"), np.empty(0, "f8")
        return np.concatenate(dts), np.concatenate(ts)

    @cached_property
    def domain_ind(self):
        return self.oct_handler.domain_ind(self.selector)

    def deposit(self, positions, fields=None, method=None, kernel_name="cubic"):
        r"""Operate on the mesh, in a particle-against-mesh fashion, with
        exclusively local input.

        This uses the octree indexing system to call a "deposition" operation
        (defined in yt/geometry/particle_deposit.pyx) that can take input from
        several particles (local to the mesh) and construct some value on the
        mesh.  The canonical example is to sum the total mass in a mesh cell
        and then divide by its volume.

        Parameters
        ----------
        positions : array_like (Nx3)
            The positions of all of the particles to be examined.  A new
            indexed octree will be constructed on these particles.
        fields : list of arrays
            All the necessary fields for computing the particle operation.  For
            instance, this might include mass, velocity, etc.
        method : string
            This is the "method name" which will be looked up in the
            `particle_deposit` namespace as `methodname_deposit`.  Current
            methods include `count`, `simple_smooth`, `sum`, `std`, `cic`,
            `weighted_mean`, `mesh_id`, and `nearest`.
        kernel_name : string, default 'cubic'
            This is the name of the smoothing kernel to use. Current supported
            kernel names include `cubic`, `quartic`, `quintic`, `wendland2`,
            `wendland4`, and `wendland6`.

        Returns
        -------
        List of fortran-ordered, mesh-like arrays.
        """
        # Here we perform our particle deposition.
        if fields is None:
            fields = []
        cls = getattr(particle_deposit, f"deposit_{method}", None)
        if cls is None:
            raise YTParticleDepositionNotImplemented(method)
        nz = self.nz
        nvals = (nz, nz, nz, (self.domain_ind >= 0).sum())
        if np.max(self.domain_ind) >= nvals[-1]:
            print(
                "nocts, domain_ind >= 0, max {} {} {}".format(
                    self.oct_handler.nocts, nvals[-1], np.max(self.domain_ind)
                )
            )
            raise Exception()
        # We allocate number of zones, not number of octs
        op = cls(nvals, kernel_name)
        op.initialize()
        mylog.debug(
            "Depositing %s (%s^3) particles into %s Octs",
            positions.shape[0],
            positions.shape[0] ** 0.3333333,
            nvals[-1],
        )
        positions.convert_to_units("code_length")
        pos = positions.d
        # We should not need the following if we know in advance all our fields
        # need no casting.
        fields = [np.ascontiguousarray(f, dtype="float64") for f in fields]
        op.process_octree(
            self.oct_handler,
            self.domain_ind,
            pos,
            fields,
            self.domain_id,
            self._domain_offset,
        )
        vals = op.finalize()
        if vals is None:
            return
        return np.asfortranarray(vals)

    def mesh_sampling_particle_field(self, positions, mesh_field, lvlmax=None):
        r"""Operate on the particles, in a mesh-against-particle
        fashion, with exclusively local input.

        This uses the octree indexing system to call a "mesh
        sampling" operation (defined in yt/geometry/particle_deposit.pyx).
        For each particle, the function returns the value of the cell
        containing the particle.

        Parameters
        ----------
        positions : array_like (Nx3)
            The positions of all of the particles to be examined.
        mesh_field : array_like (M,)
            The value of the field to deposit.
        lvlmax : array_like (N), optional
            If provided, the maximum level where to look for cells

        Returns
        -------
        List of fortran-ordered, particle-like arrays containing the
        value of the mesh at the location of the particles.
        """
        # Here we perform our particle deposition.
        npart = positions.shape[0]
        nocts = (self.domain_ind >= 0).sum()
        # We allocate number of zones, not number of octs
        op = particle_deposit.CellIdentifier(npart, "none")
        op.initialize(npart)
        mylog.debug(
            "Depositing %s Octs onto %s (%s^3) particles",
            nocts,
            positions.shape[0],
            positions.shape[0] ** 0.3333333,
        )
        pos = positions.to("code_length").value.astype("float64")

        op.process_octree(
            self.oct_handler,
            self.domain_ind,
            pos,
            None,
            self.domain_id,
            self._domain_offset,
            lvlmax=lvlmax,
        )

        igrid, icell = op.finalize()

        igrid = np.asfortranarray(igrid)
        icell = np.asfortranarray(icell)

        # Some particle may fall outside of the local domain, so we
        # need to be careful here
        ids = igrid * 8 + icell
        mask = ids >= 0

        # Fill the return array
        ret = np.empty(npart)
        ret[mask] = mesh_field[ids[mask]]
        ret[~mask] = np.nan

        return ret

    def smooth(
        self,
        positions,
        fields=None,
        index_fields=None,
        method=None,
        create_octree=False,
        nneighbors=64,
        kernel_name="cubic",
    ):
        r"""Operate on the mesh, in a particle-against-mesh fashion, with
        non-local input.

        This uses the octree indexing system to call a "smoothing" operation
        (defined in yt/geometry/particle_smooth.pyx) that can take input from
        several (non-local) particles and construct some value on the mesh.
        The canonical example is to conduct a smoothing kernel operation on the
        mesh.

        Parameters
        ----------
        positions : array_like (Nx3)
            The positions of all of the particles to be examined.  A new
            indexed octree will be constructed on these particles.
        fields : list of arrays
            All the necessary fields for computing the particle operation.  For
            instance, this might include mass, velocity, etc.
        index_fields : list of arrays
            All of the fields defined on the mesh that may be used as input to
            the operation.
        method : string
            This is the "method name" which will be looked up in the
            `particle_smooth` namespace as `methodname_smooth`.  Current
            methods include `volume_weighted`, `nearest`, `idw`,
            `nth_neighbor`, and `density`.
        create_octree : bool
            Should we construct a new octree for indexing the particles?  In
            cases where we are applying an operation on a subset of the
            particles used to construct the mesh octree, this will ensure that
            we are able to find and identify all relevant particles.
        nneighbors : int, default 64
            The number of neighbors to examine during the process.
        kernel_name : string, default 'cubic'
            This is the name of the smoothing kernel to use. Current supported
            kernel names include `cubic`, `quartic`, `quintic`, `wendland2`,
            `wendland4`, and `wendland6`.

        Returns
        -------
        List of fortran-ordered, mesh-like arrays.
        """
        # Here we perform our particle deposition.
        positions.convert_to_units("code_length")
        if create_octree:
            morton = compute_morton(
                positions[:, 0],
                positions[:, 1],
                positions[:, 2],
                self.ds.domain_left_edge,
                self.ds.domain_right_edge,
            )
            morton.sort()
            particle_octree = ParticleOctreeContainer(
                [1, 1, 1],
                self.ds.domain_left_edge,
                self.ds.domain_right_edge,
                num_zones=self._nz,
            )
            # This should ensure we get everything within one neighbor of home.
            particle_octree.n_ref = nneighbors * 2
            particle_octree.add(morton)
            particle_octree.finalize(self.domain_id)
            pdom_ind = particle_octree.domain_ind(self.selector)
        else:
            particle_octree = self.oct_handler
            pdom_ind = self.domain_ind
        if fields is None:
            fields = []
        if index_fields is None:
            index_fields = []
        cls = getattr(particle_smooth, f"{method}_smooth", None)
        if cls is None:
            raise YTParticleDepositionNotImplemented(method)
        nz = self.nz
        mdom_ind = self.domain_ind
        nvals = (nz, nz, nz, (mdom_ind >= 0).sum())
        op = cls(nvals, len(fields), nneighbors, kernel_name)
        op.initialize()
        mylog.debug(
            "Smoothing %s particles into %s Octs", positions.shape[0], nvals[-1]
        )
        # Pointer operations within 'process_octree' require arrays to be
        # contiguous cf. https://bitbucket.org/yt_analysis/yt/issues/1079
        fields = [np.ascontiguousarray(f, dtype="float64") for f in fields]
        op.process_octree(
            self.oct_handler,
            mdom_ind,
            positions,
            self.fcoords,
            fields,
            self.domain_id,
            self._domain_offset,
            self.ds.periodicity,
            index_fields,
            particle_octree,
            pdom_ind,
            self.ds.geometry,
        )
        # If there are 0s in the smoothing field this will not throw an error,
        # but silently return nans for vals where dividing by 0
        # Same as what is currently occurring, but suppressing the div by zero
        # error.
        with np.errstate(invalid="ignore"):
            vals = op.finalize()
        if isinstance(vals, list):
            vals = [np.asfortranarray(v) for v in vals]
        else:
            vals = np.asfortranarray(vals)
        return vals

    def particle_operation(
        self, positions, fields=None, method=None, nneighbors=64, kernel_name="cubic"
    ):
        r"""Operate on particles, in a particle-against-particle fashion.

        This uses the octree indexing system to call a "smoothing" operation
        (defined in yt/geometry/particle_smooth.pyx) that expects to be called
        in a particle-by-particle fashion.  For instance, the canonical example
        of this would be to compute the Nth nearest neighbor, or to compute the
        density for a given particle based on some kernel operation.

        Many of the arguments to this are identical to those used in the smooth
        and deposit functions.  Note that the `fields` argument must not be
        empty, as these fields will be modified in place.

        Parameters
        ----------
        positions : array_like (Nx3)
            The positions of all of the particles to be examined.  A new
            indexed octree will be constructed on these particles.
        fields : list of arrays
            All the necessary fields for computing the particle operation.  For
            instance, this might include mass, velocity, etc.  One of these
            will likely be modified in place.
        method : string
            This is the "method name" which will be looked up in the
            `particle_smooth` namespace as `methodname_smooth`.
        nneighbors : int, default 64
            The number of neighbors to examine during the process.
        kernel_name : string, default 'cubic'
            This is the name of the smoothing kernel to use. Current supported
            kernel names include `cubic`, `quartic`, `quintic`, `wendland2`,
            `wendland4`, and `wendland6`.

        Returns
        -------
        Nothing.

        """
        # Here we perform our particle deposition.
        positions.convert_to_units("code_length")
        morton = compute_morton(
            positions[:, 0],
            positions[:, 1],
            positions[:, 2],
            self.ds.domain_left_edge,
            self.ds.domain_right_edge,
        )
        morton.sort()
        particle_octree = ParticleOctreeContainer(
            [1, 1, 1],
            self.ds.domain_left_edge,
            self.ds.domain_right_edge,
            num_zones=2,
        )
        particle_octree.n_ref = nneighbors * 2
        particle_octree.add(morton)
        particle_octree.finalize()
        pdom_ind = particle_octree.domain_ind(self.selector)
        if fields is None:
            fields = []
        cls = getattr(particle_smooth, f"{method}_smooth", None)
        if cls is None:
            raise YTParticleDepositionNotImplemented(method)
        nz = self.nz
        mdom_ind = self.domain_ind
        nvals = (nz, nz, nz, (mdom_ind >= 0).sum())
        op = cls(nvals, len(fields), nneighbors, kernel_name)
        op.initialize()
        mylog.debug(
            "Smoothing %s particles into %s Octs", positions.shape[0], nvals[-1]
        )
        op.process_particles(
            particle_octree,
            pdom_ind,
            positions,
            fields,
            self.domain_id,
            self._domain_offset,
            self.ds.periodicity,
            self.ds.geometry,
        )
        vals = op.finalize()
        if vals is None:
            return
        if isinstance(vals, list):
            vals = [np.asfortranarray(v) for v in vals]
        else:
            vals = np.asfortranarray(vals)
        return vals

    @cell_count_cache
    def select_icoords(self, dobj):
        return self.oct_handler.icoords(
            dobj.selector, domain_id=self.domain_id, num_cells=self._cell_count
        )

    @cell_count_cache
    def select_fcoords(self, dobj):
        fcoords = self.oct_handler.fcoords(
            dobj.selector, domain_id=self.domain_id, num_cells=self._cell_count
        )
        return self.ds.arr(fcoords, "code_length")

    @cell_count_cache
    def select_fwidth(self, dobj):
        fwidth = self.oct_handler.fwidth(
            dobj.selector, domain_id=self.domain_id, num_cells=self._cell_count
        )
        return self.ds.arr(fwidth, "code_length")

    @cell_count_cache
    def select_ires(self, dobj):
        return self.oct_handler.ires(
            dobj.selector, domain_id=self.domain_id, num_cells=self._cell_count
        )

    def select(self, selector, source, dest, offset):
        n = self.oct_handler.selector_fill(
            selector, source, dest, offset, domain_id=self.domain_id
        )
        return n

    def count(self, selector):
        return -1

    def count_particles(self, selector, x, y, z):
        # We don't cache the selector results
        count = selector.count_points(x, y, z, 0.0)
        return count

    def select_particles(self, selector, x, y, z):
        mask = selector.select_points(x, y, z, 0.0)
        return mask

    def get_vertex_centered_data(self, fields):
        # Make sure the field list has only unique entries
        fields = list(set(fields))
        new_fields = {}
        cg = self.retrieve_ghost_zones(1, fields)
        for field in fields:
            new_fields[field] = cg[field][1:, 1:, 1:].copy()
            np.add(new_fields[field], cg[field][:-1, 1:, 1:], new_fields[field])
            np.add(new_fields[field], cg[field][1:, :-1, 1:], new_fields[field])
            np.add(new_fields[field], cg[field][1:, 1:, :-1], new_fields[field])
            np.add(new_fields[field], cg[field][:-1, 1:, :-1], new_fields[field])
            np.add(new_fields[field], cg[field][1:, :-1, :-1], new_fields[field])
            np.add(new_fields[field], cg[field][:-1, :-1, 1:], new_fields[field])
            np.add(new_fields[field], cg[field][:-1, :-1, :-1], new_fields[field])
            np.multiply(new_fields[field], 0.125, new_fields[field])

        return new_fields


class OctreeSubsetBlockSlicePosition:
    def __init__(self, ind, block_slice):
        self.ind = ind
        self.block_slice = block_slice
        nz = self.block_slice.octree_subset.nz
        self.ActiveDimensions = np.array([nz, nz, nz], dtype="int64")
        self.ds = block_slice.ds

    def __getitem__(self, key):
        bs = self.block_slice
        rv = np.require(
            bs.octree_subset[key][:, :, :, self.ind].T,
            requirements=bs.octree_subset._block_order,
        )
        return rv

    @property
    def id(self):
        return self.ind

    @property
    def Level(self):
        return self.block_slice._ires[0, 0, 0, self.ind]

    @property
    def LeftEdge(self):
        LE = (
            self.block_slice._fcoords[0, 0, 0, self.ind, :].d
            - self.block_slice._fwidth[0, 0, 0, self.ind, :].d * 0.5
        )
        return self.block_slice.octree_subset.ds.arr(
            LE, self.block_slice._fcoords.units
        )

    @property
    def RightEdge(self):
        RE = (
            self.block_slice._fcoords[-1, -1, -1, self.ind, :].d
            + self.block_slice._fwidth[-1, -1, -1, self.ind, :].d * 0.5
        )
        return self.block_slice.octree_subset.ds.arr(
            RE, self.block_slice._fcoords.units
        )

    @property
    def dds(self):
        return self.block_slice._fwidth[0, 0, 0, self.ind, :]

    def clear_data(self):
        pass

    def get_vertex_centered_data(self, fields, smoothed=False, no_ghost=False):
        field = fields[0]
        new_field = self.block_slice.get_vertex_centered_data(fields)[field]
        return {field: new_field[..., self.ind]}

    @contextmanager
    def _field_parameter_state(self, field_parameters):
        yield self.block_slice.octree_subset._field_parameter_state(field_parameters)


class OctreeSubsetBlockSlice:
    def __init__(self, octree_subset, ds):
        self.octree_subset = octree_subset
        self.ds = ds
        self._vertex_centered_data = {}
        # Cache some attributes
        for attr in ["ires", "icoords", "fcoords", "fwidth"]:
            v = getattr(octree_subset, attr)
            setattr(self, f"_{attr}", octree_subset._reshape_vals(v))

    @property
    def octree_subset_with_gz(self):
        subset_with_gz = getattr(self, "_octree_subset_with_gz", None)
        if not subset_with_gz:
            self._octree_subset_with_gz = self.octree_subset.retrieve_ghost_zones(1, [])
        return self._octree_subset_with_gz

    def get_vertex_centered_data(self, fields, smoothed=False, no_ghost=False):
        if no_ghost is True:
            raise NotImplementedError(
                "get_vertex_centered_data without ghost zones for "
                "oct-based datasets has not been implemented."
            )

        # Make sure the field list has only unique entries
        fields = list(set(fields))
        new_fields = {}
        cg = self.octree_subset_with_gz
        for field in fields:
            if field in self._vertex_centered_data:
                new_fields[field] = self._vertex_centered_data[field]
            else:
                finfo = self.ds._get_field_info(field)
                orig_field = cg[field]
                nocts = orig_field.shape[-1]
                new_field = np.zeros((3, 3, 3, nocts), order="F")

                # Compute vertex-centred data as mean of 8 neighbours cell data
                slices = (slice(1, None), slice(None, -1))
                for slx, sly, slz in product(*repeat(slices, 3)):
                    new_field += orig_field[slx, sly, slz]
                new_field *= 0.125

                new_fields[field] = self.ds.arr(new_field, finfo.output_units)

                self._vertex_centered_data[field] = new_fields[field]

        return new_fields

    def __iter__(self):
        for i in range(self._ires.shape[-1]):
            yield i, OctreeSubsetBlockSlicePosition(i, self)


class YTPositionArray(unyt_array):
    @property
    def morton(self):
        self.validate()
        eps = np.finfo(self.dtype).eps
        LE = self.min(axis=0)
        LE -= np.abs(LE) * eps
        RE = self.max(axis=0)
        RE += np.abs(RE) * eps
        morton = compute_morton(self[:, 0], self[:, 1], self[:, 2], LE, RE)
        return morton

    def to_octree(self, num_zones=2, dims=(1, 1, 1), n_ref=64):
        mi = self.morton
        mi.sort()
        eps = np.finfo(self.dtype).eps
        LE = self.min(axis=0)
        LE -= np.abs(LE) * eps
        RE = self.max(axis=0)
        RE += np.abs(RE) * eps
        octree = ParticleOctreeContainer(dims, LE, RE, num_zones=num_zones)
        octree.n_ref = n_ref
        octree.add(mi)
        octree.finalize()
        return octree

    def validate(self):
        if (
            len(self.shape) != 2
            or self.shape[1] != 3
            or self.units.dimensions != length
        ):
            raise YTInvalidPositionArray(self.shape, self.units.dimensions)
