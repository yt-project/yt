"""
Subsets of octrees




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.data_objects.data_containers import \
    YTFieldData, \
    YTSelectionContainer
import yt.geometry.particle_deposit as particle_deposit
import yt.geometry.particle_smooth as particle_smooth

from yt.funcs import mylog
from yt.utilities.lib.geometry_utils import compute_morton
from yt.geometry.particle_oct_container import \
    ParticleOctreeContainer
from yt.units.yt_array import YTArray
from yt.units.dimensions import length
from yt.utilities.exceptions import \
    YTInvalidPositionArray, \
    YTFieldTypeNotFound, \
    YTParticleDepositionNotImplemented

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
    _type_name = 'octree_subset'
    _skip_add = True
    _con_args = ('base_region', 'domain', 'ds')
    _domain_offset = 0
    _cell_count = -1
    _block_reorder = None

    def __init__(self, base_region, domain, ds, over_refine_factor = 1):
        self._num_zones = 1 << (over_refine_factor)
        self._oref = over_refine_factor
        self.field_data = YTFieldData()
        self.field_parameters = {}
        self.domain = domain
        self.domain_id = domain.domain_id
        self.ds = domain.ds
        self._index = self.ds.index
        self.oct_handler = domain.oct_handler
        self._last_mask = None
        self._last_selector_id = None
        self._current_particle_type = 'all'
        self._current_fluid_type = self.ds.default_fluid_type
        self.base_region = base_region
        self.base_selector = base_region.selector

    def __getitem__(self, key):
        tr = super(OctreeSubset, self).__getitem__(key)
        try:
            fields = self._determine_fields(key)
        except YTFieldTypeNotFound:
            return tr
        finfo = self.ds._get_field_info(*fields[0])
        if not finfo.particle_type:
            # We may need to reshape the field, if it is being queried from
            # field_data.  If it's already cached, it just passes through.
            if len(tr.shape) < 4:
                tr = self._reshape_vals(tr)
            return tr
        return tr

    @property
    def nz(self):
        return self._num_zones + 2*self._num_ghost_zones

    def _reshape_vals(self, arr):
        nz = self.nz
        if len(arr.shape) <= 2:
            n_oct = arr.shape[0] // (nz**3)
        else:
            n_oct = max(arr.shape)
        if arr.size == nz*nz*nz*n_oct:
            new_shape = (nz, nz, nz, n_oct)
        elif arr.size == nz*nz*nz*n_oct * 3:
            new_shape = (nz, nz, nz, n_oct, 3)
        else:
            raise RuntimeError
        # Note that if arr is already F-contiguous, this *shouldn't* copy the
        # data.  But, it might.  However, I can't seem to figure out how to
        # make the assignment to .shape, which *won't* copy the data, make the
        # resultant array viewed in Fortran order.
        arr = arr.reshape(new_shape, order="F")
        return arr

    _domain_ind = None

    def mask_refinement(self, selector):
        mask = self.oct_handler.mask(selector, domain_id = self.domain_id)
        return mask

    def select_blocks(self, selector):
        mask = self.oct_handler.mask(selector, domain_id = self.domain_id)
        slicer = OctreeSubsetBlockSlice(self)
        for i, sl in slicer:
            yield sl, np.atleast_3d(mask[i,...])

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

    @property
    def domain_ind(self):
        if self._domain_ind is None:
            di = self.oct_handler.domain_ind(self.selector)
            self._domain_ind = di
        return self._domain_ind

    def deposit(self, positions, fields = None, method = None,
                kernel_name='cubic'):
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
        if fields is None: fields = []
        cls = getattr(particle_deposit, "deposit_%s" % method, None)
        if cls is None:
            raise YTParticleDepositionNotImplemented(method)
        nz = self.nz
        nvals = (nz, nz, nz, (self.domain_ind >= 0).sum())
        # We allocate number of zones, not number of octs
        op = cls(nvals, kernel_name)
        op.initialize()
        mylog.debug("Depositing %s (%s^3) particles into %s Octs",
            positions.shape[0], positions.shape[0]**0.3333333, nvals[-1])
        pos = np.asarray(positions.convert_to_units("code_length"),
                         dtype="float64")
        # We should not need the following if we know in advance all our fields
        # need no casting.
        fields = [np.asarray(f, dtype="float64") for f in fields]
        op.process_octree(self.oct_handler, self.domain_ind, pos, fields,
            self.domain_id, self._domain_offset)
        vals = op.finalize()
        if vals is None: return
        return np.asfortranarray(vals)

    def smooth(self, positions, fields = None, index_fields = None,
               method = None, create_octree = False, nneighbors = 64,
               kernel_name = 'cubic'):
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
                positions[:,0], positions[:,1], positions[:,2],
                self.ds.domain_left_edge,
                self.ds.domain_right_edge)
            morton.sort()
            particle_octree = ParticleOctreeContainer([1, 1, 1],
                self.ds.domain_left_edge,
                self.ds.domain_right_edge,
                over_refine = self._oref)
            # This should ensure we get everything within one neighbor of home.
            particle_octree.n_ref = nneighbors * 2
            particle_octree.add(morton)
            particle_octree.finalize()
            pdom_ind = particle_octree.domain_ind(self.selector)
        else:
            particle_octree = self.oct_handler
            pdom_ind = self.domain_ind
        if fields is None: fields = []
        if index_fields is None: index_fields = []
        cls = getattr(particle_smooth, "%s_smooth" % method, None)
        if cls is None:
            raise YTParticleDepositionNotImplemented(method)
        nz = self.nz
        mdom_ind = self.domain_ind
        nvals = (nz, nz, nz, (mdom_ind >= 0).sum())
        op = cls(nvals, len(fields), nneighbors, kernel_name)
        op.initialize()
        mylog.debug("Smoothing %s particles into %s Octs",
            positions.shape[0], nvals[-1])
        # Pointer operations within 'process_octree' require arrays to be
        # contiguous cf. https://bitbucket.org/yt_analysis/yt/issues/1079
        fields = [np.ascontiguousarray(f, dtype="float64") for f in fields]
        op.process_octree(self.oct_handler, mdom_ind, positions,
            self.fcoords, fields,
            self.domain_id, self._domain_offset, self.ds.periodicity,
            index_fields, particle_octree, pdom_ind, self.ds.geometry)
        # If there are 0s in the smoothing field this will not throw an error,
        # but silently return nans for vals where dividing by 0
        # Same as what is currently occurring, but suppressing the div by zero
        # error.
        with np.errstate(invalid='ignore'):
            vals = op.finalize()
        if vals is None: return
        if isinstance(vals, list):
            vals = [np.asfortranarray(v) for v in vals]
        else:
            vals = np.asfortranarray(vals)
        return vals

    def particle_operation(self, positions, fields = None,
            method = None, nneighbors = 64, kernel_name = 'cubic'):
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
            positions[:,0], positions[:,1], positions[:,2],
            self.ds.domain_left_edge,
            self.ds.domain_right_edge)
        morton.sort()
        particle_octree = ParticleOctreeContainer([1, 1, 1],
            self.ds.domain_left_edge,
            self.ds.domain_right_edge,
            over_refine = 1)
        particle_octree.n_ref = nneighbors * 2
        particle_octree.add(morton)
        particle_octree.finalize()
        pdom_ind = particle_octree.domain_ind(self.selector)
        if fields is None: fields = []
        cls = getattr(particle_smooth, "%s_smooth" % method, None)
        if cls is None:
            raise YTParticleDepositionNotImplemented(method)
        nz = self.nz
        mdom_ind = self.domain_ind
        nvals = (nz, nz, nz, (mdom_ind >= 0).sum())
        op = cls(nvals, len(fields), nneighbors, kernel_name)
        op.initialize()
        mylog.debug("Smoothing %s particles into %s Octs",
            positions.shape[0], nvals[-1])
        op.process_particles(particle_octree, pdom_ind, positions,
            fields, self.domain_id, self._domain_offset, self.ds.periodicity,
            self.ds.geometry)
        vals = op.finalize()
        if vals is None: return
        if isinstance(vals, list):
            vals = [np.asfortranarray(v) for v in vals]
        else:
            vals = np.asfortranarray(vals)
        return vals

    @cell_count_cache
    def select_icoords(self, dobj):
        return self.oct_handler.icoords(dobj.selector, domain_id = self.domain_id,
                                     num_cells = self._cell_count)

    @cell_count_cache
    def select_fcoords(self, dobj):
        fcoords = self.oct_handler.fcoords(
            dobj.selector, domain_id=self.domain_id, num_cells=self._cell_count)
        return self.ds.arr(fcoords, 'code_length')

    @cell_count_cache
    def select_fwidth(self, dobj):
        fwidth = self.oct_handler.fwidth(
            dobj.selector, domain_id=self.domain_id, num_cells=self._cell_count)
        return self.ds.arr(fwidth, 'code_length')

    @cell_count_cache
    def select_ires(self, dobj):
        return self.oct_handler.ires(dobj.selector, domain_id = self.domain_id,
                                     num_cells = self._cell_count)

    def select(self, selector, source, dest, offset):
        n = self.oct_handler.selector_fill(selector, source, dest, offset,
                                           domain_id = self.domain_id)
        return n

    def count(self, selector):
        return -1

    def count_particles(self, selector, x, y, z):
        # We don't cache the selector results
        count = selector.count_points(x,y,z, 0.0)
        return count

    def select_particles(self, selector, x, y, z):
        mask = selector.select_points(x,y,z, 0.0)
        return mask

class ParticleOctreeSubset(OctreeSubset):
    # Subclassing OctreeSubset is somewhat dubious.
    # This is some subset of an octree.  Note that the sum of subsets of an
    # octree may multiply include data files.  While we can attempt to mitigate
    # this, it's unavoidable for many types of data storage on disk.
    _type_name = 'indexed_octree_subset'
    _con_args = ('data_files', 'ds', 'min_ind', 'max_ind')
    domain_id = -1
    def __init__(self, base_region, data_files, ds, min_ind = 0, max_ind = 0,
                 over_refine_factor = 1):
        # The first attempt at this will not work in parallel.
        self._num_zones = 1 << (over_refine_factor)
        self._oref = over_refine_factor
        self.data_files = data_files
        self.field_data = YTFieldData()
        self.field_parameters = {}
        self.ds = ds
        self._index = self.ds.index
        self.oct_handler = ds.index.oct_handler
        self.min_ind = min_ind
        if max_ind == 0: max_ind = (1 << 63)
        self.max_ind = max_ind
        self._last_mask = None
        self._last_selector_id = None
        self._current_particle_type = 'all'
        self._current_fluid_type = self.ds.default_fluid_type
        self.base_region = base_region
        self.base_selector = base_region.selector

class OctreeSubsetBlockSlicePosition(object):
    def __init__(self, ind, block_slice):
        self.ind = ind
        self.block_slice = block_slice
        nz = self.block_slice.octree_subset.nz
        self.ActiveDimensions = np.array([nz,nz,nz], dtype="int64")

    def __getitem__(self, key):
        bs = self.block_slice
        rv = bs.octree_subset[key][:,:,:,self.ind]
        if bs.octree_subset._block_reorder:
            rv = rv.copy(order=bs.octree_subset._block_reorder)
        return rv

    @property
    def id(self):
        return self.ind

    @property
    def Level(self):
        return self.block_slice._ires[0,0,0,self.ind]

    @property
    def LeftEdge(self):
        LE = (self.block_slice._fcoords[0,0,0,self.ind,:]
            - self.block_slice._fwidth[0,0,0,self.ind,:]*0.5)
        return LE

    @property
    def RightEdge(self):
        RE = (self.block_slice._fcoords[-1,-1,-1,self.ind,:]
            + self.block_slice._fwidth[-1,-1,-1,self.ind,:]*0.5)
        return RE

    @property
    def dds(self):
        return self.block_slice._fwidth[0,0,0,self.ind,:]

    def clear_data(self):
        pass

    def get_vertex_centered_data(self, *args, **kwargs):
        raise NotImplementedError


class OctreeSubsetBlockSlice(object):
    def __init__(self, octree_subset):
        self.octree_subset = octree_subset
        # Cache some attributes
        for attr in ["ires", "icoords", "fcoords", "fwidth"]:
            v = getattr(octree_subset, attr)
            setattr(self, "_%s" % attr, octree_subset._reshape_vals(v))

    def __iter__(self):
        for i in range(self._ires.shape[-1]):
            yield i, OctreeSubsetBlockSlicePosition(i, self)

class YTPositionArray(YTArray):
    @property
    def morton(self):
        self.validate()
        eps = np.finfo(self.dtype).eps
        LE = self.min(axis=0)
        LE -= np.abs(LE) * eps
        RE = self.max(axis=0)
        RE += np.abs(RE) * eps
        morton = compute_morton(
            self[:,0], self[:,1], self[:,2],
            LE, RE)
        return morton

    def to_octree(self, over_refine_factor = 1, dims = (1,1,1),
                  n_ref = 64):
        mi = self.morton
        mi.sort()
        eps = np.finfo(self.dtype).eps
        LE = self.min(axis=0)
        LE -= np.abs(LE) * eps
        RE = self.max(axis=0)
        RE += np.abs(RE) * eps
        octree = ParticleOctreeContainer(dims, LE, RE,
            over_refine = over_refine_factor)
        octree.n_ref = n_ref
        octree.add(mi)
        octree.finalize()
        return octree

    def validate(self):
        if len(self.shape) != 2 or self.shape[1] != 3 \
           or self.units.dimensions != length:
            raise YTInvalidPositionArray(self.shape, self.units.dimensions)
