# distutils: include_dirs = LIB_DIR_EWAH
# distutils: language = c++
# distutils: extra_compile_args = CPP14_FLAG
# distutils: libraries = STD_LIBS
"""
Oct container tuned for Particles




"""


from libc.math cimport ceil, floor, fmod, log2
from libc.stdlib cimport free, malloc, qsort
from libc.string cimport memset
from libcpp.map cimport map as cmap
from libcpp.vector cimport vector

from yt.utilities.lib.ewah_bool_array cimport (
    bool_array,
    ewah_bool_array,
    ewah_bool_iterator,
    ewah_map,
    ewah_word_type,
)

import numpy as np

cimport cython
cimport numpy as np
cimport oct_visitors
from cpython.exc cimport PyErr_CheckSignals
from cython cimport floating
from cython.operator cimport dereference, preincrement
from oct_container cimport (
    ORDER_MAX,
    Oct,
    OctAllocationContainer,
    OctInfo,
    OctKey,
    OctreeContainer,
    SparseOctreeContainer,
)
from oct_visitors cimport OctVisitor, cind
from selection_routines cimport AlwaysSelector, SelectorObject

from yt.utilities.lib.fp_utils cimport *
from yt.utilities.lib.geometry_utils cimport (
    bounded_morton,
    bounded_morton_dds,
    bounded_morton_relative_dds,
    bounded_morton_split_dds,
    bounded_morton_split_relative_dds,
    decode_morton_64bit,
    encode_morton_64bit,
    morton_neighbors_coarse,
    morton_neighbors_refined,
)

from collections import defaultdict

from yt.funcs import get_pbar

from particle_deposit cimport gind

#from yt.utilities.lib.ewah_bool_wrap cimport \
from ..utilities.lib.ewah_bool_wrap cimport BoolArrayCollection

import os
import struct

# If set to 1, ghost cells are added at the refined level reguardless of if the
# coarse cell containing it is refined in the selector.
# If set to 0, ghost cells are only added at the refined level if the coarse
# index for the ghost cell is refined in the selector.
DEF RefinedExternalGhosts = 1

_bitmask_version = np.uint64(5)

from ..utilities.lib.ewah_bool_wrap cimport (
    BoolArrayCollectionUncompressed as BoolArrayColl,
    FileBitmasks,
    SparseUnorderedBitmaskSet as SparseUnorderedBitmask,
    SparseUnorderedRefinedBitmaskSet as SparseUnorderedRefinedBitmask,
)

ctypedef cmap[np.uint64_t, bool_array] CoarseRefinedSets

cdef class ParticleOctreeContainer(OctreeContainer):
    cdef Oct** oct_list
    #The starting oct index of each domain
    cdef np.int64_t *dom_offsets
    cdef public int max_level
    #How many particles do we keep befor refining
    cdef public int n_ref

    def allocate_root(self):
        cdef int i, j, k
        cdef Oct *cur
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                for k in range(self.nn[2]):
                    cur = self.allocate_oct()
                    self.root_mesh[i][j][k] = cur

    def __dealloc__(self):
        #Call the freemem ops on every ocy
        #of the root mesh recursively
        cdef int i, j, k
        if self.root_mesh == NULL: return
        for i in range(self.nn[0]):
            if self.root_mesh[i] == NULL: continue
            for j in range(self.nn[1]):
                if self.root_mesh[i][j] == NULL: continue
                for k in range(self.nn[2]):
                    if self.root_mesh[i][j][k] == NULL: continue
                    self.visit_free(self.root_mesh[i][j][k])
        free(self.oct_list)
        free(self.dom_offsets)

    cdef void visit_free(self, Oct *o):
        #Free the memory for this oct recursively
        cdef int i, j, k
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if o.children != NULL \
                       and o.children[cind(i,j,k)] != NULL:
                        self.visit_free(o.children[cind(i,j,k)])
        free(o.children)
        free(o)

    def clear_fileind(self):
        cdef int i, j, k
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                for k in range(self.nn[2]):
                    self.visit_clear(self.root_mesh[i][j][k])

    cdef void visit_clear(self, Oct *o):
        #Free the memory for this oct recursively
        cdef int i, j, k
        o.file_ind = 0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if o.children != NULL \
                       and o.children[cind(i,j,k)] != NULL:
                        self.visit_clear(o.children[cind(i,j,k)])

    def __iter__(self):
        #Get the next oct, will traverse domains
        #Note that oct containers can be sorted
        #so that consecutive octs are on the same domain
        cdef int oi
        cdef Oct *o
        for oi in range(self.nocts):
            o = self.oct_list[oi]
            yield (o.file_ind, o.domain_ind, o.domain)

    def allocate_domains(self, domain_counts):
        pass

    def finalize(self, int domain_id = 0):
        #This will sort the octs in the oct list
        #so that domains appear consecutively
        #And then find the oct index/offset for
        #every domain
        cdef int max_level = 0
        self.oct_list = <Oct**> malloc(sizeof(Oct*)*self.nocts)
        cdef np.int64_t i = 0, lpos = 0
        # Note that we now assign them in the same order they will be visited
        # by recursive visitors.
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                for k in range(self.nn[2]):
                    self.visit_assign(self.root_mesh[i][j][k], &lpos,
                                      0, &max_level)
        assert(lpos == self.nocts)
        for i in range(self.nocts):
            self.oct_list[i].domain_ind = i
            self.oct_list[i].domain = domain_id
        self.max_level = max_level

    cdef visit_assign(self, Oct *o, np.int64_t *lpos, int level, int *max_level):
        cdef int i, j, k
        self.oct_list[lpos[0]] = o
        lpos[0] += 1
        max_level[0] = imax(max_level[0], level)
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if o.children != NULL \
                       and o.children[cind(i,j,k)] != NULL:
                        self.visit_assign(o.children[cind(i,j,k)], lpos,
                                level + 1, max_level)
        return

    cdef np.int64_t get_domain_offset(self, int domain_id):
        return 0

    cdef Oct* allocate_oct(self):
        #Allocate the memory, set to NULL or -1
        #We reserve space for n_ref particles, but keep
        #track of how many are used with np initially 0
        self.nocts += 1
        cdef Oct *my_oct = <Oct*> malloc(sizeof(Oct))
        my_oct.domain = -1
        my_oct.file_ind = 0
        my_oct.domain_ind = self.nocts - 1
        my_oct.children = NULL
        return my_oct

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def add(self, np.ndarray[np.uint64_t, ndim=1] indices,
            np.uint8_t order = ORDER_MAX):
        #Add this particle to the root oct
        #Then if that oct has children, add it to them recursively
        #If the child needs to be refined because of max particles, do so
        cdef np.int64_t no = indices.shape[0], p
        cdef np.uint64_t index
        cdef int i, level
        cdef int ind[3]
        if self.root_mesh[0][0][0] == NULL: self.allocate_root()
        cdef np.uint64_t *data = <np.uint64_t *> indices.data
        for p in range(no):
            # We have morton indices, which means we choose left and right by
            # looking at (MAX_ORDER - level) & with the values 1, 2, 4.
            level = 0
            index = indices[p]
            if index == FLAG:
                # This is a marker for the index not being inside the domain
                # we're interested in.
                continue
            # Convert morton index to 3D index of octree root
            for i in range(3):
                ind[i] = (index >> ((order - level)*3 + (2 - i))) & 1
            cur = self.root_mesh[ind[0]][ind[1]][ind[2]]
            if cur == NULL:
                raise RuntimeError
            # Continue refining the octree until you reach the level of the
            # morton indexing order. Along the way, use prefix to count
            # previous indices at levels in the octree?
            while (cur.file_ind + 1) > self.n_ref:
                if level >= order: break # Just dump it here.
                level += 1
                for i in range(3):
                    ind[i] = (index >> ((order - level)*3 + (2 - i))) & 1
                if cur.children == NULL or \
                   cur.children[cind(ind[0],ind[1],ind[2])] == NULL:
                    cur = self.refine_oct(cur, index, level, order)
                    self.filter_particles(cur, data, p, level, order)
                else:
                    cur = cur.children[cind(ind[0],ind[1],ind[2])]
            # If our n_ref is 1, we are always refining, which means we're an
            # index octree.  In this case, we should store the index for fast
            # lookup later on when we find neighbors and the like.
            if self.n_ref == 1:
                cur.file_ind = index
            else:
                cur.file_ind += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Oct *refine_oct(self, Oct *o, np.uint64_t index, int level,
                         np.uint8_t order):
        #Allocate and initialize child octs
        #Attach particles to child octs
        #Remove particles from this oct entirely
        cdef int i, j, k
        cdef int ind[3]
        cdef Oct *noct
        # TODO: This does not need to be changed.
        o.children = <Oct **> malloc(sizeof(Oct *)*8)
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    noct = self.allocate_oct()
                    noct.domain = o.domain
                    noct.file_ind = 0
                    o.children[cind(i,j,k)] = noct
        o.file_ind = self.n_ref + 1
        for i in range(3):
            ind[i] = (index >> ((order - level)*3 + (2 - i))) & 1
        noct = o.children[cind(ind[0],ind[1],ind[2])]
        return noct

    cdef void filter_particles(self, Oct *o, np.uint64_t *data, np.int64_t p,
                               int level, np.uint8_t order):
        # Now we look at the last nref particles to decide where they go.
        # If p: Loops over all previous morton indices
        # If n_ref: Loops over n_ref previous morton indices
        cdef int n = imin(p, self.n_ref)
        cdef np.uint64_t *arr = data + imax(p - self.n_ref, 0)
        cdef np.uint64_t prefix1, prefix2
        # Now we figure out our prefix, which is the oct address at this level.
        # As long as we're actually in Morton order, we do not need to worry
        # about *any* of the other children of the oct.
        prefix1 = data[p] >> (order - level)*3
        for i in range(n):
            prefix2 = arr[i] >> (order - level)*3
            if (prefix1 == prefix2):
                o.file_ind += 1 # Says how many morton indices are in this octant?

    def recursively_count(self):
        #Visit every cell, accumulate the # of cells per level
        cdef int i, j, k
        cdef np.int64_t counts[128]
        for i in range(128): counts[i] = 0
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                for k in range(self.nn[2]):
                    if self.root_mesh[i][j][k] != NULL:
                        self.visit(self.root_mesh[i][j][k], counts)
        level_counts = {}
        for i in range(128):
            if counts[i] == 0: break
            level_counts[i] = counts[i]
        return level_counts

    cdef visit(self, Oct *o, np.int64_t *counts, level = 0):
        cdef int i, j, k
        counts[level] += 1
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if o.children != NULL \
                       and o.children[cind(i,j,k)] != NULL:
                        self.visit(o.children[cind(i,j,k)], counts, level + 1)
        return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Oct *get_from_index(self, np.uint64_t mi, np.uint8_t order = ORDER_MAX,
                             int max_level = 99):
        cdef Oct *cur
        cdef Oct *next
        cur = next = NULL
        cdef int i
        cdef np.int64_t level = -1
        cdef int ind32[3]
        cdef np.uint64_t ind[3]
        cdef np.uint64_t index
        # Get level offset
        cdef int level_offset[3]
        for i in range(3):
            level_offset[i] = np.log2(self.nn[i])
            if (1 << level_offset[i]) != self.nn[i]:
                raise Exception("Octree does not have roots along dimension {} in a power of 2 ".format(i))
        for i in range(2,3):
            if level_offset[i] != level_offset[0]:
                raise Exception("Octree must have the same number of roots in each dimension for this.")
        # Get root for index
        index = (mi >> ((order - level_offset[0])*3))
        decode_morton_64bit(index, ind)
        for i in range(3):
            ind32[i] = ind[i]
        self.get_root(ind32, &next)
        # We want to stop recursing when there's nowhere else to go
        level = level_offset[0]
        max_level = min(max_level, order)
        while next != NULL and level <= max_level:
            level += 1
            for i in range(3):
                ind[i] = (mi >> ((order - level)*3 + (2 - i))) & 1
            cur = next
            if cur.children != NULL:
                next = cur.children[cind(ind[0],ind[1],ind[2])]
            else:
                next = NULL
        return cur

    def apply_domain(self, int domain_id, BoolArrayCollection mask,
                     int masklevel):
        cdef SelectorObject selector = AlwaysSelector(None)
        ind = self.domain_ind(selector, mask = mask, masklevel = masklevel)
        for i in range(self.nocts):
            if ind[i] < 0: continue
            self.oct_list[i].domain = domain_id
        ind_out = super(ParticleOctreeContainer,self).domain_ind(selector, domain_id = domain_id)

    def domain_ind(self, selector, int domain_id = -1,
                   BoolArrayCollection mask = None, int masklevel = 99):
        if mask is None:
            return super(ParticleOctreeContainer,self).domain_ind(selector, domain_id = domain_id)
        # Create mask for octs that are touched by the mask
        cdef ewah_bool_array *ewah_slct = <ewah_bool_array *> mask.ewah_keys
        cdef ewah_bool_iterator *iter_set = new ewah_bool_iterator(ewah_slct[0].begin())
        cdef ewah_bool_iterator *iter_end = new ewah_bool_iterator(ewah_slct[0].end())
        cdef np.ndarray[np.uint8_t, ndim=1] oct_mask
        oct_mask = np.zeros(self.nocts, 'uint8')
        cdef Oct *o
        cdef int coct, cmi
        coct = cmi = 0
        while iter_set[0] != iter_end[0]:
            mi = dereference(iter_set[0])
            o = self.get_from_index(mi, order = masklevel)
            if o != NULL:
                _mask_children(oct_mask, o)
                coct += 1
            cmi += 1;
            preincrement(iter_set[0])
        # Get domain ind
        cdef np.ndarray[np.int64_t, ndim=1] ind
        ind = np.zeros(self.nocts, 'int64') - 1
        cdef oct_visitors.MaskedIndexOcts visitor
        visitor = oct_visitors.MaskedIndexOcts(self, domain_id)
        visitor.oct_index = ind
        visitor.oct_mask = oct_mask
        self.visit_all_octs(selector, visitor)
        return ind

cdef void _mask_children(np.ndarray[np.uint8_t] mask, Oct *cur):
    cdef int i, j, k
    if cur == NULL:
        return
    mask[cur.domain_ind] = 1
    if cur.children == NULL:
        return
    for i in range(2):
        for j in range(2):
            for k in range(2):
                _mask_children(mask, cur.children[cind(i,j,k)])

cdef np.uint64_t ONEBIT=1
cdef np.uint64_t FLAG = ~(<np.uint64_t>0)

cdef class ParticleBitmap:
    cdef np.float64_t left_edge[3]
    cdef np.float64_t right_edge[3]
    cdef np.uint8_t periodicity[3]
    cdef np.float64_t dds[3]
    cdef np.float64_t dds_mi1[3]
    cdef np.float64_t dds_mi2[3]
    cdef np.float64_t idds[3]
    cdef np.int32_t dims[3]
    cdef np.int64_t file_hash
    cdef np.uint64_t directional_max2[3]
    cdef public np.uint64_t nfiles
    cdef public np.int32_t index_order1
    cdef public np.int32_t index_order2
    cdef public object masks
    cdef public object particle_counts
    cdef public object counts
    cdef public object max_count
    cdef public object _last_selector
    cdef public object _last_return_values
    cdef public object _cached_octrees
    cdef public object _last_octree_subset
    cdef public object _last_oct_handler
    cdef public object _prev_octree_subset
    cdef public object _prev_oct_handler
    cdef np.uint32_t *file_markers
    cdef np.uint64_t n_file_markers
    cdef np.uint64_t file_marker_i
    cdef public FileBitmasks bitmasks
    cdef public BoolArrayCollection collisions
    cdef public int _used_mi2

    def __init__(self, left_edge, right_edge, periodicity, file_hash, nfiles,
                 index_order1, index_order2):
        # TODO: Set limit on maximum orders?
        cdef int i
        self._cached_octrees = {}
        self._last_selector = None
        self._last_return_values = None
        self._last_octree_subset = None
        self._last_oct_handler = None
        self._prev_octree_subset = None
        self._prev_oct_handler = None
        self.file_hash = file_hash
        self.nfiles = nfiles
        for i in range(3):
            self.left_edge[i] = left_edge[i]
            self.right_edge[i] = right_edge[i]
            self.periodicity[i] = <np.uint8_t>periodicity[i]
            self.dims[i] = (1<<index_order1)
            self.dds[i] = (right_edge[i] - left_edge[i])/self.dims[i]
            self.idds[i] = 1.0/self.dds[i]
            self.dds_mi1[i] = (right_edge[i] - left_edge[i]) / (1<<index_order1)
        # We use 64-bit masks
        self._used_mi2 = 0
        self.index_order1 = index_order1
        self._update_mi2(index_order2)
        # This will be an on/off flag for which morton index values are touched
        # by particles.
        # This is the simple way, for now.
        self.masks = np.zeros((1 << (index_order1 * 3), nfiles), dtype="uint8")
        self.particle_counts = np.zeros(1 << (index_order1 * 3), dtype="uint64")
        self.bitmasks = FileBitmasks(self.nfiles)
        self.collisions = BoolArrayCollection()

    def _bitmask_logicaland(self, ifile, bcoll, out):
        self.bitmasks._logicaland(ifile, bcoll, out)

    def _bitmask_intersects(self, ifile, bcoll):
        return self.bitmasks._intersects(ifile, bcoll)

    def update_mi2(self, np.float64_t characteristic_size,
                   np.uint64_t max_index_order2 = 6):
        """
        mi2 is the *refined* morton index order; mi2 is thus the definition of
        the size of the refined index objects we stick inside any collisions at
        the coarse level.  This takes a characteristic size and attempts to
        compute the mi2 such that the cell is roughly equivalent to the
        characteristic size.  It will return whether or not it was able to
        update; if the mi2 has already been used, it does not update.
        There are cases where the maximum index_order2 that it would compute
        would be extremely fine, which can do really bad things to memory.  So
        we allow the setting of a maximum value, which it won't exceed.
        """
        if self._used_mi2 > 0:
            return self.index_order2
        cdef np.uint64_t index_order2 = 2
        for i in range(3):
            # Note we're casting to signed here, to avoid negative issues.
            if self.dds_mi1[i] < characteristic_size: continue
            index_order2 = max(index_order2, <np.uint64_t> ceil(log2(self.dds_mi1[i] / characteristic_size)))
        index_order2 = i64min(max_index_order2, index_order2)
        self._update_mi2(index_order2)
        return self.index_order2

    cdef void _update_mi2(self, np.uint64_t index_order2):
        self.index_order2 = index_order2
        mi2_max = (1 << self.index_order2) - 1
        self.directional_max2[0] = encode_morton_64bit(mi2_max, 0, 0)
        self.directional_max2[1] = encode_morton_64bit(0, mi2_max, 0)
        self.directional_max2[2] = encode_morton_64bit(0, 0, mi2_max)
        for i in range(3):
            self.dds_mi2[i] = self.dds_mi1[i] / (1<<index_order2)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def _coarse_index_data_file(self,
                                np.ndarray[floating, ndim=2] pos,
                                np.ndarray[floating, ndim=1] hsml,
                                np.uint64_t file_id):
        return self.__coarse_index_data_file(pos, hsml, file_id)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void __coarse_index_data_file(self,
                                       np.ndarray[floating, ndim=2] pos,
                                       np.ndarray[floating, ndim=1] hsml,
                                       np.uint64_t file_id) except *:
        # Initialize
        cdef np.int64_t i, p
        cdef np.uint64_t mi, miex
        cdef np.uint64_t mi_split[3]
        cdef np.float64_t ppos[3]
        cdef np.float64_t s_ppos[3] # shifted ppos
        cdef np.float64_t clip_pos_l[3]
        cdef np.float64_t clip_pos_r[3]
        cdef int skip
        cdef np.uint64_t bounds[2][3]
        cdef np.uint64_t xex, yex, zex
        cdef np.float64_t LE[3]
        cdef np.float64_t RE[3]
        cdef np.float64_t DW[3]
        cdef np.uint8_t PER[3]
        cdef np.float64_t dds[3]
        cdef np.float64_t radius
        cdef np.uint8_t[:] mask = self.masks[:, file_id]
        cdef np.uint64_t[:] particle_counts = self.particle_counts
        cdef np.uint64_t msize = (1 << (self.index_order1 * 3))
        cdef int axiter[3][2]
        cdef np.float64_t axiterv[3][2]
        # Copy over things for this file (type cast necessary?)
        for i in range(3):
            LE[i] = self.left_edge[i]
            RE[i] = self.right_edge[i]
            PER[i] = self.periodicity[i]
            dds[i] = self.dds_mi1[i]
            DW[i] = RE[i] - LE[i]
            axiter[i][0] = 0 # We always do an offset of 0
            axiterv[i][0] = 0.0
        # Mark index of particles that are in this file
        for p in range(pos.shape[0]):
            skip = 0
            for i in range(3):
                axiter[i][1] = 999
                # Skip particles outside the domain
                if not (LE[i] <= pos[p, i] < RE[i]):
                    skip = 1
                    break
                ppos[i] = pos[p,i]
            if skip == 1: continue
            mi = bounded_morton_split_dds(ppos[0], ppos[1], ppos[2], LE,
                                          dds, mi_split)
            mask[mi] = 1
            particle_counts[mi] += 1
            # Expand mask by softening
            if hsml is None:
                continue
            if hsml[p] < 0:
                raise RuntimeError(
                    f"Smoothing length for particle {p} is negative with "
                    f"value {hsml[p]}")
            radius = hsml[p]
            # We first check if we're bounded within the domain; this follows the logic in the
            # pixelize_cartesian routine.  We assume that no smoothing
            # length can wrap around both directions.
            for i in range(3):
                if PER[i] and ppos[i] - radius < LE[i]:
                    axiter[i][1] = +1
                    axiterv[i][1] = DW[i]
                elif PER[i] and ppos[i] + radius > RE[i]:
                    axiter[i][1] = -1
                    axiterv[i][1] = -DW[i]
            for xi in range(2):
                if axiter[0][xi] == 999: continue
                s_ppos[0] = ppos[0] + axiterv[0][xi]
                for yi in range(2):
                    if axiter[1][yi] == 999: continue
                    s_ppos[1] = ppos[1] + axiterv[1][yi]
                    for zi in range(2):
                        if axiter[2][zi] == 999: continue
                        s_ppos[2] = ppos[2] + axiterv[2][zi]
                        # OK, now we compute the left and right edges for this shift.
                        for i in range(3):
                            clip_pos_l[i] = fmax(s_ppos[i] - radius, LE[i] + dds[i]/10)
                            clip_pos_r[i] = fmin(s_ppos[i] + radius, RE[i] - dds[i]/10)
                        bounded_morton_split_dds(clip_pos_l[0], clip_pos_l[1], clip_pos_l[2], LE, dds, bounds[0])
                        bounded_morton_split_dds(clip_pos_r[0], clip_pos_r[1], clip_pos_r[2], LE, dds, bounds[1])
                        # We go to the upper bound plus one so that we have *inclusive* loops -- the upper bound
                        # is the cell *index*, so we want to make sure we include that cell.  This is also why
                        # we don't need to worry about mi_max being the max index rather than the cell count.
                        for xex in range(bounds[0][0], bounds[1][0] + 1):
                            for yex in range(bounds[0][1], bounds[1][1] + 1):
                                for zex in range(bounds[0][2], bounds[1][2] + 1):
                                    miex = encode_morton_64bit(xex, yex, zex)
                                    mask[miex] = 1
                                    particle_counts[miex] += 1
                                    if miex >= msize:
                                        raise IndexError(
                                            "Index for a softening region "
                                            f"({miex}) exceeds "
                                            f"max ({msize})")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def _set_coarse_index_data_file(self, np.uint64_t file_id):
        return self.__set_coarse_index_data_file(file_id)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void __set_coarse_index_data_file(self, np.uint64_t file_id):
        cdef np.int64_t i
        cdef FileBitmasks bitmasks = self.bitmasks
        cdef np.ndarray[np.uint8_t, ndim=1] mask = self.masks[:,file_id]
        # Add in order
        for i in range(mask.shape[0]):
            if mask[i] == 1:
                bitmasks._set_coarse(file_id, i)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    def _refined_index_data_file(self,
                                 BoolArrayCollection in_collection,
                                 np.ndarray[floating, ndim=2] pos,
                                 np.ndarray[floating, ndim=1] hsml,
                                 np.ndarray[np.uint8_t, ndim=1] mask,
                                 np.ndarray[np.uint64_t, ndim=1] sub_mi1,
                                 np.ndarray[np.uint64_t, ndim=1] sub_mi2,
                                 np.uint64_t file_id, np.int64_t nsub_mi,
                                 np.uint64_t count_threshold = 128,
                                 np.uint8_t mask_threshold = 2):
        self._used_mi2 = 1
        if in_collection is None:
            in_collection = BoolArrayCollection()
        cdef BoolArrayCollection _in_coll = in_collection
        out_collection = self.__refined_index_data_file(_in_coll, pos, hsml, mask,
                                              count_threshold, mask_threshold)
        return 0, out_collection

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef BoolArrayCollection __refined_index_data_file(
        self,
        BoolArrayCollection in_collection,
        np.ndarray[floating, ndim=2] pos,
        np.ndarray[floating, ndim=1] hsml,
        np.ndarray[np.uint8_t, ndim=1] mask,
        np.uint64_t count_threshold, np.uint8_t mask_threshold
    ):
        # Initialize
        cdef np.int64_t p, sorted_ind
        cdef np.uint64_t i
        cdef np.uint64_t mi1, mi2
        cdef np.float64_t ppos[3]
        cdef np.float64_t s_ppos[3] # shifted ppos
        cdef int skip
        cdef BoolArrayCollection this_collection, out_collection
        cdef np.uint64_t bounds[2][3]
        cdef np.uint8_t fully_enclosed
        cdef np.float64_t LE[3]
        cdef np.float64_t RE[3]
        cdef np.float64_t DW[3]
        cdef np.uint8_t PER[3]
        cdef np.float64_t dds1[3]
        cdef np.float64_t dds2[3]
        cdef np.float64_t radius
        cdef np.uint64_t mi_split1[3]
        cdef np.uint64_t mi_split2[3]
        cdef np.uint64_t miex1
        cdef np.uint64_t[:] particle_counts = self.particle_counts
        cdef np.uint64_t xex, yex, zex
        cdef np.float64_t clip_pos_l[3]
        cdef np.float64_t clip_pos_r[3]
        cdef int axiter[3][2]
        cdef np.float64_t axiterv[3][2]
        cdef CoarseRefinedSets coarse_refined_map
        cdef cmap[np.uint64_t, np.uint64_t] refined_count
        cdef np.uint64_t nfully_enclosed = 0, n_calls = 0
        mi1_max = (1 << self.index_order1) - 1
        mi2_max = (1 << self.index_order2) - 1
        cdef np.uint64_t max_mi1_elements = 1 << (3*self.index_order1)
        cdef np.uint64_t max_mi2_elements = 1 << (3*self.index_order2)
        for i in range(max_mi1_elements):
            refined_count[i] = 0
        # Copy things from structure (type cast)
        for i in range(3):
            LE[i] = self.left_edge[i]
            RE[i] = self.right_edge[i]
            PER[i] = self.periodicity[i]
            dds1[i] = self.dds_mi1[i]
            dds2[i] = self.dds_mi2[i]
            DW[i] = RE[i] - LE[i]
            axiter[i][0] = 0 # We always do an offset of 0
            axiterv[i][0] = 0.0
        cdef np.ndarray[np.uint64_t, ndim=1] morton_indices = np.empty(pos.shape[0], dtype="u8")
        for p in range(pos.shape[0]):
            morton_indices[p] = bounded_morton(pos[p, 0], pos[p, 1], pos[p, 2],
                                               LE, RE, self.index_order1)
        # Loop over positions skipping those outside the domain
        cdef np.ndarray[np.uint64_t, ndim=1, cast=True] sorted_order
        if hsml is None:
            # casting to uint64 for compatibility with 32 bits systems
            # see https://github.com/yt-project/yt/issues/3656
            sorted_order = np.argsort(morton_indices).astype(np.uint64, copy=False)
        else:
            sorted_order = np.argsort(hsml)[::-1].astype(np.uint64, copy=False)
        for sorted_ind in range(sorted_order.shape[0]):
            p = sorted_order[sorted_ind]
            skip = 0
            for i in range(3):
                axiter[i][1] = 999
                if not (LE[i] <= pos[p, i] < RE[i]):
                    skip = 1
                    break
                ppos[i] = pos[p,i]
            if skip == 1: continue
            # Only look if collision at coarse index
            mi1 = bounded_morton_split_dds(ppos[0], ppos[1], ppos[2], LE,
                                           dds1, mi_split1)
            if hsml is None:
                if mask[mi1] < mask_threshold \
                        or particle_counts[mi1] < count_threshold:
                    continue
                # Determine sub index within cell of primary index
                mi2 = bounded_morton_split_relative_dds(
                    ppos[0], ppos[1], ppos[2], LE, dds1, dds2, mi_split2)
                if refined_count[mi1] == 0:
                    coarse_refined_map[mi1].padWithZeroes(max_mi2_elements)
                if not coarse_refined_map[mi1].get(mi2):
                    coarse_refined_map[mi1].set(mi2)
                    refined_count[mi1] += 1
            else: # only hit if we have smoothing lengths.
                # We have to do essentially the identical process to in the coarse indexing,
                # except here we need to fill in all the subranges as well as the coarse ranges
                # Note that we are also doing the null case, where we do no shifting
                radius = hsml[p]
                #if mask[mi1] <= 4: # only one thing in this area
                #    continue
                for i in range(3):
                    if PER[i] and ppos[i] - radius < LE[i]:
                        axiter[i][1] = +1
                        axiterv[i][1] = DW[i]
                    elif PER[i] and ppos[i] + radius > RE[i]:
                        axiter[i][1] = -1
                        axiterv[i][1] = -DW[i]
                for xi in range(2):
                    if axiter[0][xi] == 999: continue
                    s_ppos[0] = ppos[0] + axiterv[0][xi]
                    for yi in range(2):
                        if axiter[1][yi] == 999: continue
                        s_ppos[1] = ppos[1] + axiterv[1][yi]
                        for zi in range(2):
                            if axiter[2][zi] == 999: continue
                            s_ppos[2] = ppos[2] + axiterv[2][zi]
                            # OK, now we compute the left and right edges for this shift.
                            for i in range(3):
                                # casting to int64 is not nice but is so we can have negative values we clip
                                clip_pos_l[i] = fmax(s_ppos[i] - radius, LE[i] + dds1[i]/10)
                                clip_pos_r[i] = fmin(s_ppos[i] + radius, RE[i] - dds1[i]/10)

                            bounded_morton_split_dds(clip_pos_l[0], clip_pos_l[1], clip_pos_l[2], LE, dds1, bounds[0])
                            bounded_morton_split_dds(clip_pos_r[0], clip_pos_r[1], clip_pos_r[2], LE, dds1, bounds[1])

                            # We go to the upper bound plus one so that we have *inclusive* loops -- the upper bound
                            # is the cell *index*, so we want to make sure we include that cell.  This is also why
                            # we don't need to worry about mi_max being the max index rather than the cell count.
                            # One additional thing to note is that for all of
                            # the *internal* cells, i.e., those that are both
                            # greater than the left edge and less than the
                            # right edge, we are fully enclosed.
                            for xex in range(bounds[0][0], bounds[1][0] + 1):
                                for yex in range(bounds[0][1], bounds[1][1] + 1):
                                    for zex in range(bounds[0][2], bounds[1][2] + 1):
                                        miex1 = encode_morton_64bit(xex, yex, zex)
                                        if mask[miex1] < mask_threshold or \
                                                particle_counts[miex1] < count_threshold:
                                            continue
                                        # this explicitly requires that it be *between*
                                        # them, not overlapping
                                        if xex > bounds[0][0] and xex < bounds[1][0] and \
                                           yex > bounds[0][1] and yex < bounds[1][1] and \
                                           zex > bounds[0][2] and zex < bounds[1][2]:
                                            fully_enclosed = 1
                                        else:
                                            fully_enclosed = 0
                                        # Now we need to fill our sub-range
                                        if refined_count[miex1] == 0:
                                            coarse_refined_map[miex1].padWithZeroes(max_mi2_elements)
                                        elif refined_count[miex1] >= max_mi2_elements:
                                            continue
                                        if fully_enclosed == 1:
                                            nfully_enclosed += 1
                                            coarse_refined_map[miex1].inplace_logicalxor(
                                                coarse_refined_map[miex1])
                                            coarse_refined_map[miex1].inplace_logicalnot()
                                            refined_count[miex1] = max_mi2_elements
                                            continue
                                        n_calls += 1
                                        refined_count[miex1] += self.__fill_refined_ranges(s_ppos, radius, LE, RE,
                                                                   dds1, xex, yex, zex,
                                                                   dds2,
                                                                   coarse_refined_map[miex1])
        cdef np.uint64_t vec_i
        cdef bool_array *buf = NULL
        cdef ewah_word_type w
        this_collection = BoolArrayCollection()
        cdef ewah_bool_array *refined_arr = NULL
        for it1 in coarse_refined_map:
            mi1 = it1.first
            refined_arr = &this_collection.ewah_coll[0][mi1]
            this_collection.ewah_keys[0].set(mi1)
            this_collection.ewah_refn[0].set(mi1)
            buf = &it1.second
            for vec_i in range(buf.sizeInBytes() / sizeof(ewah_word_type)):
                w = buf.getWord(vec_i)
                refined_arr.addWord(w)
        out_collection = BoolArrayCollection()
        in_collection._logicalor(this_collection, out_collection)
        return out_collection

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef np.int64_t __fill_refined_ranges(self, np.float64_t s_ppos[3], np.float64_t radius,
                                           np.float64_t LE[3], np.float64_t RE[3],
                                           np.float64_t dds1[3], np.uint64_t xex, np.uint64_t yex, np.uint64_t zex,
                                           np.float64_t dds2[3], bool_array &refined_set) except -1:
        cdef int i
        cdef np.uint64_t bounds_l[3], bounds_r[3]
        cdef np.uint64_t miex2, miex2_min, miex2_max
        cdef np.float64_t clip_pos_l[3]
        cdef np.float64_t clip_pos_r[3]
        cdef np.float64_t cell_edge_l, cell_edge_r
        cdef np.uint64_t ex1[3]
        cdef np.uint64_t xiex_min, yiex_min, ziex_min
        cdef np.uint64_t xiex_max, yiex_max, ziex_max
        cdef np.uint64_t old_nsub = refined_set.numberOfOnes()
        ex1[0] = xex
        ex1[1] = yex
        ex1[2] = zex
        # Check a few special cases
        for i in range(3):
            # Figure out our bounds inside our coarse cell, in the space of the
            # full domain
            cell_edge_l = ex1[i] * dds1[i] + LE[i]
            cell_edge_r = cell_edge_l + dds1[i]
            if s_ppos[i] + radius < cell_edge_l or s_ppos[i] - radius > cell_edge_r:
                return 0
            clip_pos_l[i] = fmax(s_ppos[i] - radius, cell_edge_l + dds2[i]/2.0)
            clip_pos_r[i] = fmin(s_ppos[i] + radius, cell_edge_r - dds2[i]/2.0)
        miex2_min = bounded_morton_split_relative_dds(clip_pos_l[0], clip_pos_l[1], clip_pos_l[2],
                                                LE, dds1, dds2, bounds_l)
        miex2_max = bounded_morton_split_relative_dds(clip_pos_r[0], clip_pos_r[1], clip_pos_r[2],
                                                LE, dds1, dds2, bounds_r)
        xex_max = self.directional_max2[0]
        yex_max = self.directional_max2[1]
        zex_max = self.directional_max2[2]
        xiex_min = miex2_min & xex_max
        yiex_min = miex2_min & yex_max
        ziex_min = miex2_min & zex_max
        xiex_max = miex2_max & xex_max
        yiex_max = miex2_max & yex_max
        ziex_max = miex2_max & zex_max
        # This could *probably* be sped up by iterating over words.
        for miex2 in range(miex2_min, miex2_max + 1):
            #miex2 = encode_morton_64bit(xex2, yex2, zex2)
            #decode_morton_64bit(miex2, ex2)
            # Let's check all our cases here
            if (miex2 & xex_max) < (xiex_min): continue
            if (miex2 & xex_max) > (xiex_max): continue
            if (miex2 & yex_max) < (yiex_min): continue
            if (miex2 & yex_max) > (yiex_max): continue
            if (miex2 & zex_max) < (ziex_min): continue
            if (miex2 & zex_max) > (ziex_max): continue
            refined_set.set(miex2)
        return refined_set.numberOfOnes() - old_nsub

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    def _set_refined_index_data_file(self,
                                     np.ndarray[np.uint64_t, ndim=1] sub_mi1,
                                     np.ndarray[np.uint64_t, ndim=1] sub_mi2,
                                     np.uint64_t file_id, np.int64_t nsub_mi):
        return self.__set_refined_index_data_file(sub_mi1, sub_mi2,
                                                  file_id, nsub_mi)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void __set_refined_index_data_file(self,
                                            np.ndarray[np.uint64_t, ndim=1] sub_mi1,
                                            np.ndarray[np.uint64_t, ndim=1] sub_mi2,
                                            np.uint64_t file_id, np.int64_t nsub_mi):
        cdef np.int64_t i, p
        cdef FileBitmasks bitmasks = self.bitmasks
        bitmasks._set_refined_index_array(file_id, nsub_mi, sub_mi1, sub_mi2)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    def find_collisions(self, verbose=False):
        cdef tuple cc, rc
        cc, rc = self.bitmasks._find_collisions(self.collisions,verbose)
        return cc, rc

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    def find_collisions_coarse(self, verbose=False, file_list = None):
        cdef int nc, nm
        nc, nm = self.bitmasks._find_collisions_coarse(self.collisions, verbose, file_list)
        return nc, nm

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    def find_uncontaminated(self, np.uint32_t ifile, BoolArrayCollection mask,
                            BoolArrayCollection mask2 = None):
        cdef np.ndarray[np.uint8_t, ndim=1] arr = np.zeros((1 << (self.index_order1 * 3)),'uint8')
        cdef np.uint8_t[:] arr_view = arr
        self.bitmasks._select_uncontaminated(ifile, mask, arr_view, mask2)
        return arr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    def find_contaminated(self, np.uint32_t ifile, BoolArrayCollection mask,
                          BoolArrayCollection mask2 = None):
        cdef np.ndarray[np.uint8_t, ndim=1] arr = np.zeros((1 << (self.index_order1 * 3)),'uint8')
        cdef np.uint8_t[:] arr_view = arr
        cdef np.ndarray[np.uint8_t, ndim=1] sfiles = np.zeros(self.nfiles,'uint8')
        cdef np.uint8_t[:] sfiles_view = sfiles
        self.bitmasks._select_contaminated(ifile, mask, arr_view, sfiles_view, mask2)
        return arr, np.where(sfiles)[0].astype('uint32')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    def find_collisions_refined(self, verbose=False):
        cdef np.int32_t nc, nm
        nc, nm = self.bitmasks._find_collisions_refined(self.collisions,verbose)
        return nc, nm

    def calcsize_bitmasks(self):
        # TODO: All cython
        cdef bytes serial_BAC
        cdef np.uint64_t ifile
        cdef int out = 0
        out += struct.calcsize('Q')
        # Bitmaps for each file
        for ifile in range(self.nfiles):
            serial_BAC = self.bitmasks._dumps(ifile)
            out += struct.calcsize('Q')
            out += len(serial_BAC)
        # Bitmap for collisions
        serial_BAC = self.collisions._dumps()
        out += struct.calcsize('Q')
        out += len(serial_BAC)
        return out

    def get_bitmasks(self):
        return self.bitmasks

    def iseq_bitmask(self, solf):
        return self.bitmasks._iseq(solf.get_bitmasks())

    def save_bitmasks(self, fname):
        cdef bytes serial_BAC
        cdef np.uint64_t ifile
        f = open(fname,'wb')
        # Header
        f.write(struct.pack('Q', _bitmask_version))
        f.write(struct.pack('q', self.file_hash))
        f.write(struct.pack('Q', self.nfiles))
        # Bitmap for each file
        for ifile in range(self.nfiles):
            serial_BAC = self.bitmasks._dumps(ifile)
            f.write(struct.pack('Q', len(serial_BAC)))
            f.write(serial_BAC)
        # Collisions
        serial_BAC = self.collisions._dumps()
        f.write(struct.pack('Q', len(serial_BAC)))
        f.write(serial_BAC)
        f.close()

    def check_bitmasks(self):
        return self.bitmasks._check()

    def reset_bitmasks(self):
        self.bitmasks._reset()

    def load_bitmasks(self, fname):
        cdef bint read_flag = 1
        cdef bint irflag
        cdef np.uint64_t ver
        cdef np.uint64_t nfiles = 0
        cdef np.int64_t file_hash
        cdef np.uint64_t size_serial
        cdef bint overwrite = 0
        # Verify that file is correct version
        if not os.path.isfile(fname):
            raise OSError
        f = open(fname,'rb')
        ver, = struct.unpack('Q',f.read(struct.calcsize('Q')))
        if ver == self.nfiles and ver != _bitmask_version:
            overwrite = 1
            nfiles = ver
            ver = 0 # Original bitmaps had number of files first
        if ver != _bitmask_version:
            raise OSError("The file format of the index has changed since "
                          "this file was created. It will be replaced with an "
                          "updated version.")
        # Read file hash
        file_hash, = struct.unpack('q', f.read(struct.calcsize('q')))
        if file_hash != self.file_hash:
            raise OSError
        # Read number of bitmaps
        if nfiles == 0:
            nfiles, = struct.unpack('Q', f.read(struct.calcsize('Q')))
            if nfiles != self.nfiles:
                raise OSError(
                    "Number of bitmasks ({}) conflicts with number of files "
                    "({})".format(nfiles, self.nfiles))
        # Read bitmap for each file
        pb = get_pbar("Loading particle index", nfiles)
        for ifile in range(nfiles):
            pb.update(ifile+1)
            size_serial, = struct.unpack('Q', f.read(struct.calcsize('Q')))
            irflag = self.bitmasks._loads(ifile, f.read(size_serial))
            if irflag == 0:
                read_flag = 0
        pb.finish()
        # Collisions
        size_serial, = struct.unpack('Q',f.read(struct.calcsize('Q')))
        irflag = self.collisions._loads(f.read(size_serial))
        f.close()
        # Save in correct format
        if overwrite == 1:
            self.save_bitmasks(fname)
        return read_flag

    def print_info(self):
        cdef np.uint64_t ifile
        for ifile in range(self.nfiles):
            self.bitmasks.print_info(ifile, "File: %03d" % ifile)

    def count_coarse(self, ifile):
        r"""Get the number of coarse cells set for a file."""
        return self.bitmasks.count_coarse(ifile)

    def count_refined(self, ifile):
        r"""Get the number of cells refined for a file."""
        return self.bitmasks.count_refined(ifile)

    def count_total(self, ifile):
        r"""Get the total number of cells set for a file."""
        return self.bitmasks.count_total(ifile)

    def check(self):
        cdef np.uint64_t mi1
        cdef ewah_bool_array arr_totref, arr_tottwo
        cdef ewah_bool_array arr, arr_any, arr_two, arr_swap
        cdef vector[size_t] vec_totref
        cdef vector[size_t].iterator it_mi1
        cdef int nm = 0, nc = 0
        cdef np.uint64_t ifile, nbitmasks
        nbitmasks = len(self.bitmasks)
        # Locate all indices with second level refinement
        for ifile in range(self.nfiles):
            arr = (<ewah_bool_array**> self.bitmasks.ewah_refn)[ifile][0]
            arr_totref.logicalor(arr,arr_totref)
        # Count collections & second level indices
        vec_totref = arr_totref.toArray()
        it_mi1 = vec_totref.begin()
        while it_mi1 != vec_totref.end():
            mi1 = dereference(it_mi1)
            arr_any.reset()
            arr_two.reset()
            for ifile in range(nbitmasks):
                if self.bitmasks._isref(ifile, mi1) == 1:
                    arr = (<cmap[np.int64_t, ewah_bool_array]**> self.bitmasks.ewah_coll)[ifile][0][mi1]
                    arr_any.logicaland(arr, arr_two) # Indices in previous files
                    arr_any.logicalor(arr, arr_swap) # All second level indices
                    arr_any = arr_swap
                    arr_two.logicalor(arr_tottwo,arr_tottwo)
            nc += arr_tottwo.numberOfOnes()
            nm += arr_any.numberOfOnes()
            preincrement(it_mi1)
        # nc: total number of second level morton indices that are repeated
        # nm: total number of second level morton indices
        print("Total of %s / %s collisions (% 3.5f%%)" % (nc, nm, 100.0*float(nc)/nm))

    def primary_indices(self):
        mi = (<ewah_bool_array*> self.collisions.ewah_keys)[0].toArray()
        return np.array(mi,'uint64')

    def file_ownership_mask(self, fid):
        cdef BoolArrayCollection out
        out = self.bitmasks._get_bitmask(<np.uint32_t> fid)
        return out

    def finalize(self):
        return
        # self.index_octree = ParticleOctreeContainer([1,1,1],
        #     [self.left_edge[0], self.left_edge[1], self.left_edge[2]],
        #     [self.right_edge[0], self.right_edge[1], self.right_edge[2]],
        #     over_refine = 0
        # )
        # self.index_octree.n_ref = 1
        # mi = (<ewah_bool_array*> self.collisions.ewah_keys)[0].toArray()
        # Change from vector to numpy
        # mi = mi.astype("uint64")
        # self.index_octree.add(mi, self.index_order1)
        # self.index_octree.finalize()

    def get_DLE(self):
        cdef int i
        cdef np.ndarray[np.float64_t, ndim=1] DLE
        DLE = np.zeros(3, dtype='float64')
        for i in range(3):
            DLE[i] = self.left_edge[i]
        return DLE
    def get_DRE(self):
        cdef int i
        cdef np.ndarray[np.float64_t, ndim=1] DRE
        DRE = np.zeros(3, dtype='float64')
        for i in range(3):
            DRE[i] = self.right_edge[i]
        return DRE

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def get_ghost_zones(self, SelectorObject selector, int ngz,
                        BoolArrayCollection dmask = None, bint coarse_ghosts = False):
        cdef BoolArrayCollection gmask, gmask2, out
        cdef np.ndarray[np.uint8_t, ndim=1] periodic = selector.get_periodicity()
        cdef bint periodicity[3]
        cdef int i
        for i in range(3):
            periodicity[i] = periodic[i]
        if dmask is None:
            dmask = BoolArrayCollection()
            gmask2 = BoolArrayCollection()
            morton_selector = ParticleBitmapSelector(selector,self,ngz=0)
            morton_selector.fill_masks(dmask, gmask2)
        gmask = BoolArrayCollection()
        dmask._get_ghost_zones(ngz, self.index_order1, self.index_order2,
                               periodicity, gmask, <bint>coarse_ghosts)
        dfiles, gfiles = self.masks_to_files(dmask, gmask)
        out = BoolArrayCollection()
        gmask._logicalor(dmask, out)
        return gfiles, out

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def selector2mask(self, SelectorObject selector):
        cdef BoolArrayCollection cmask = BoolArrayCollection()
        cdef ParticleBitmapSelector morton_selector
        morton_selector = ParticleBitmapSelector(selector,self,ngz=0)
        morton_selector.fill_masks(cmask)
        return cmask

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def mask2files(self, BoolArrayCollection cmask):
        cdef np.ndarray[np.uint32_t, ndim=1] file_idx
        file_idx = self.mask_to_files(cmask)
        return file_idx

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def mask2filemasks(self, BoolArrayCollection cmask, np.ndarray[np.uint32_t, ndim=1] file_idx):
        cdef BoolArrayCollection fmask
        cdef np.int32_t fid
        cdef np.ndarray[object, ndim=1] file_masks
        cdef int i
        # Get bitmasks for parts of files touching the selector
        file_masks = np.array([BoolArrayCollection() for i in range(len(file_idx))],
                              dtype="object")
        for i, (fid, fmask) in enumerate(zip(file_idx,file_masks)):
            self.bitmasks._logicaland(<np.uint32_t> fid, cmask, fmask)
        return file_masks

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def filemasks2addfiles(self, np.ndarray[object, ndim=1] file_masks):
        cdef list addfile_idx
        addfile_idx = len(file_masks)*[None]
        for i, fmask in enumerate(file_masks):
            addfile_idx[i] = self.mask_to_files(fmask).astype('uint32')
        return addfile_idx

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def identify_file_masks(self, SelectorObject selector):
        cdef BoolArrayCollection cmask = BoolArrayCollection()
        cdef BoolArrayCollection fmask
        cdef np.int32_t fid
        cdef np.ndarray[object, ndim=1] file_masks
        cdef np.ndarray[np.uint32_t, ndim=1] file_idx
        cdef list addfile_idx
        # Get bitmask for selector
        cdef ParticleBitmapSelector morton_selector
        morton_selector = ParticleBitmapSelector(selector, self, ngz=0)
        morton_selector.fill_masks(cmask)
        # Get bitmasks for parts of files touching the selector
        file_idx = self.mask_to_files(cmask)
        file_masks = np.array([BoolArrayCollection() for i in range(len(file_idx))],
                              dtype="object")
        addfile_idx = len(file_idx)*[None]
        for i, (fid, fmask) in enumerate(zip(file_idx,file_masks)):
            self.bitmasks._logicaland(<np.uint32_t> fid, cmask, fmask)
            addfile_idx[i] = self.mask_to_files(fmask).astype('uint32')
        return file_idx.astype('uint32'), file_masks, addfile_idx

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def identify_data_files(self, SelectorObject selector, int ngz = 0):
        cdef BoolArrayCollection cmask_s = BoolArrayCollection()
        cdef BoolArrayCollection cmask_g = BoolArrayCollection()
        # Find mask of selected morton indices
        cdef ParticleBitmapSelector morton_selector
        morton_selector = ParticleBitmapSelector(selector, self, ngz=ngz)
        morton_selector.fill_masks(cmask_s, cmask_g)
        return self.masks_to_files(cmask_s, cmask_g), (cmask_s, cmask_g)

    def mask_to_files(self, BoolArrayCollection mm_s):
        cdef FileBitmasks mm_d = self.bitmasks
        cdef np.uint32_t ifile
        cdef np.ndarray[np.uint8_t, ndim=1] file_mask_p
        file_mask_p = np.zeros(self.nfiles, dtype="uint8")
        # Compare with mask of particles
        for ifile in range(self.nfiles):
            # Only continue if the file is not already selected
            if file_mask_p[ifile] == 0:
                if mm_d._intersects(ifile, mm_s):
                    file_mask_p[ifile] = 1
        cdef np.ndarray[np.int32_t, ndim=1] file_idx_p
        file_idx_p = np.where(file_mask_p)[0].astype('int32')
        return file_idx_p.astype('uint32')

    def masks_to_files(self, BoolArrayCollection mm_s, BoolArrayCollection mm_g):
        cdef FileBitmasks mm_d = self.bitmasks
        cdef np.uint32_t ifile
        cdef np.ndarray[np.uint8_t, ndim=1] file_mask_p
        cdef np.ndarray[np.uint8_t, ndim=1] file_mask_g
        file_mask_p = np.zeros(self.nfiles, dtype="uint8")
        file_mask_g = np.zeros(self.nfiles, dtype="uint8")
        # Compare with mask of particles
        for ifile in range(self.nfiles):
            # Only continue if the file is not already selected
            if file_mask_p[ifile] == 0:
                if mm_d._intersects(ifile, mm_s):
                    file_mask_p[ifile] = 1
                    file_mask_g[ifile] = 0 # No intersection
                elif mm_d._intersects(ifile, mm_g):
                    file_mask_g[ifile] = 1
        cdef np.ndarray[np.int32_t, ndim=1] file_idx_p
        cdef np.ndarray[np.int32_t, ndim=1] file_idx_g
        file_idx_p = np.where(file_mask_p)[0].astype('int32')
        file_idx_g = np.where(file_mask_g)[0].astype('int32')
        return file_idx_p.astype('uint32'), file_idx_g.astype('uint32')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def construct_octree(self, index, io_handler, data_files,
                         over_refine_factor,
                         BoolArrayCollection selector_mask,
                         BoolArrayCollection base_mask = None):
        cdef np.uint64_t total_pcount
        cdef np.uint64_t i, j, k, n
        cdef int ind[3]
        cdef np.uint64_t ind64[3]
        cdef ParticleBitmapOctreeContainer octree
        cdef np.uint64_t mi, mi_root
        cdef np.ndarray pos
        cdef np.ndarray[np.float32_t, ndim=2] pos32
        cdef np.ndarray[np.float64_t, ndim=2] pos64
        cdef np.float64_t ppos[3]
        cdef np.float64_t DLE[3]
        cdef np.float64_t DRE[3]
        cdef int bitsize = 0
        for i in range(3):
            DLE[i] = self.left_edge[i]
            DRE[i] = self.right_edge[i]
        cdef np.ndarray[np.uint64_t, ndim=1] morton_ind
        # Determine cells that need to be added to the octree
        cdef np.ndarray[np.uint32_t, ndim=1] file_idx_p
        cdef np.ndarray[np.uint32_t, ndim=1] file_idx_g
        cdef np.uint64_t nroot = selector_mask._count_total()
        # Now we can actually create a sparse octree.
        octree = ParticleBitmapOctreeContainer(
            (self.dims[0], self.dims[1], self.dims[2]),
            (self.left_edge[0], self.left_edge[1], self.left_edge[2]),
            (self.right_edge[0], self.right_edge[1], self.right_edge[2]),
            nroot, over_refine_factor)
        octree.n_ref = index.dataset.n_ref
        octree.level_offset = self.index_order1
        octree.allocate_domains()
        # Add roots based on the mask
        cdef np.uint64_t croot = 0
        cdef ewah_bool_array *ewah_slct = <ewah_bool_array *> selector_mask.ewah_keys
        cdef ewah_bool_array *ewah_base
        if base_mask is not None:
            ewah_base = <ewah_bool_array *> base_mask.ewah_keys
        else:
            ewah_base = NULL
        cdef ewah_bool_iterator *iter_set = new ewah_bool_iterator(ewah_slct[0].begin())
        cdef ewah_bool_iterator *iter_end = new ewah_bool_iterator(ewah_slct[0].end())
        cdef np.ndarray[np.uint8_t, ndim=1] slct_arr
        slct_arr = np.zeros((1 << (self.index_order1 * 3)),'uint8')
        while iter_set[0] != iter_end[0]:
            mi = dereference(iter_set[0])
            if ewah_base != NULL and ewah_base[0].get(mi) == 0:
                octree._index_base_roots[croot] = 0
                slct_arr[mi] = 2
            else:
                slct_arr[mi] = 1
            decode_morton_64bit(mi, ind64)
            for j in range(3):
                ind[j] = ind64[j]
            octree.next_root(1, ind)
            croot += 1
            preincrement(iter_set[0])
        assert(croot == nroot)
        if ewah_base != NULL:
            assert(np.sum(octree._index_base_roots) == ewah_base[0].numberOfOnes())
        # Get morton indices for all particles in this file and those
        # contaminating cells it has majority control of.
        files_touched = data_files #+ buffer_files  # datafile object from ID goes here
        total_pcount = 0
        for data_file in files_touched:
            total_pcount += sum(data_file.total_particles.values())
        morton_ind = np.empty(total_pcount, dtype='uint64')
        total_pcount = 0
        cdef np.uint64_t base_pcount = 0
        for data_file in files_touched:
            # We now get our particle positions
            for pos in io_handler._yield_coordinates(data_file):
                pos32 = pos64 = None
                bitsize = 0
                if pos.dtype == np.float32:
                    pos32 = pos
                    bitsize = 32
                    for j in range(pos.shape[0]):
                        for k in range(3):
                            ppos[k] = pos32[j,k]
                        mi = bounded_morton(ppos[0], ppos[1], ppos[2],
                                            DLE, DRE, ORDER_MAX)
                        mi_root = mi >> (3*(ORDER_MAX-self.index_order1))
                        if slct_arr[mi_root] > 0:
                            morton_ind[total_pcount] = mi
                            total_pcount += 1
                            if slct_arr[mi_root] == 1:
                                base_pcount += 1
                elif pos.dtype == np.float64:
                    pos64 = pos
                    bitsize = 64
                    for j in range(pos.shape[0]):
                        for k in range(3):
                            ppos[k] = pos64[j,k]
                        mi = bounded_morton(ppos[0], ppos[1], ppos[2],
                                            DLE, DRE, ORDER_MAX)
                        mi_root = mi >> (3*(ORDER_MAX-self.index_order1))
                        if slct_arr[mi_root] > 0:
                            morton_ind[total_pcount] = mi
                            total_pcount += 1
                            if slct_arr[mi_root] == 1:
                                base_pcount += 1
                else:
                    raise RuntimeError
        morton_ind = morton_ind[:total_pcount]
        morton_ind.sort()
        octree.add(morton_ind, self.index_order1)
        octree.finalize()
        return octree

cdef class ParticleBitmapSelector:
    cdef SelectorObject selector
    cdef ParticleBitmap bitmap
    cdef np.uint32_t ngz
    cdef np.float64_t DLE[3]
    cdef np.float64_t DRE[3]
    cdef bint periodicity[3]
    cdef np.uint32_t order1
    cdef np.uint32_t order2
    cdef np.uint64_t max_index1
    cdef np.uint64_t max_index2
    cdef np.uint64_t s1
    cdef np.uint64_t s2
    cdef void* pointers[11]
    cdef np.uint64_t[:,:] ind1_n
    cdef np.uint64_t[:,:] ind2_n
    cdef np.uint32_t[:,:] neighbors
    cdef np.uint64_t[:] neighbor_list1
    cdef np.uint64_t[:] neighbor_list2
    cdef np.uint32_t nfiles
    cdef np.uint8_t[:] file_mask_p
    cdef np.uint8_t[:] file_mask_g
    # Uncompressed boolean
    cdef np.uint8_t[:] refined_select_bool
    cdef np.uint8_t[:] refined_ghosts_bool
    cdef np.uint8_t[:] coarse_select_bool
    cdef np.uint8_t[:] coarse_ghosts_bool
    cdef SparseUnorderedRefinedBitmask refined_ghosts_list
    cdef BoolArrayColl select_ewah
    cdef BoolArrayColl ghosts_ewah

    def __cinit__(self, selector, bitmap, ngz=0):
        cdef int i
        cdef np.ndarray[np.uint8_t, ndim=1] periodicity = np.zeros(3, dtype='uint8')
        cdef np.ndarray[np.float64_t, ndim=1] DLE = np.zeros(3, dtype='float64')
        cdef np.ndarray[np.float64_t, ndim=1] DRE = np.zeros(3, dtype='float64')

        self.selector = selector
        self.bitmap = bitmap
        self.ngz = ngz
        # Things from the bitmap & selector
        periodicity = selector.get_periodicity()
        DLE = bitmap.get_DLE()
        DRE = bitmap.get_DRE()
        for i in range(3):
            self.DLE[i] = DLE[i]
            self.DRE[i] = DRE[i]
            self.periodicity[i] = periodicity[i]
        self.order1 = bitmap.index_order1
        self.order2 = bitmap.index_order2
        self.nfiles = bitmap.nfiles
        self.max_index1 = <np.uint64_t>(1 << self.order1)
        self.max_index2 = <np.uint64_t>(1 << self.order2)
        self.s1 = <np.uint64_t>(1 << (self.order1*3))
        self.s2 = <np.uint64_t>(1 << (self.order2*3))

        self.neighbors = np.zeros((2*ngz+1, 3), dtype='uint32')
        self.ind1_n = np.zeros((2*ngz+1, 3), dtype='uint64')
        self.ind2_n = np.zeros((2*ngz+1, 3), dtype='uint64')
        self.neighbor_list1 = np.zeros((2*ngz+1)**3, dtype='uint64')
        self.neighbor_list2 = np.zeros((2*ngz+1)**3, dtype='uint64')
        self.file_mask_p = np.zeros(bitmap.nfiles, dtype='uint8')
        self.file_mask_g = np.zeros(bitmap.nfiles, dtype='uint8')

        self.refined_select_bool = np.zeros(self.s2, 'uint8')
        self.refined_ghosts_bool = np.zeros(self.s2, 'uint8')
        self.coarse_select_bool = np.zeros(self.s1, 'uint8')
        self.coarse_ghosts_bool = np.zeros(self.s1, 'uint8')

        self.refined_ghosts_list = SparseUnorderedRefinedBitmask()
        self.select_ewah = BoolArrayColl(self.s1, self.s2)
        self.ghosts_ewah = BoolArrayColl(self.s1, self.s2)

    def fill_masks(self, BoolArrayCollection mm_s, BoolArrayCollection mm_g = None):
        # Normal variables
        cdef int i
        cdef np.int32_t level = 0
        cdef np.uint64_t mi1
        mi1 = ~(<np.uint64_t>0)
        cdef np.float64_t pos[3]
        cdef np.float64_t dds[3]
        cdef np.uint64_t cur_ind[3]
        for i in range(3):
            cur_ind[i] = 0
            pos[i] = self.DLE[i]
            dds[i] = self.DRE[i] - self.DLE[i]
        if mm_g is None:
            mm_g = BoolArrayCollection()
        # Uncompressed version
        cdef BoolArrayColl mm_s0
        cdef BoolArrayColl mm_g0
        mm_s0 = BoolArrayColl(self.s1, self.s2)
        mm_g0 = BoolArrayColl(self.s1, self.s2)
        # Recurse
        cdef np.float64_t rpos[3]
        for i in range(3):
            rpos[i] = self.DRE[i] - self.bitmap.dds_mi2[i]/2.0
        sbbox = self.selector.select_bbox_edge(pos, rpos)
        if sbbox == 1:
            for mi1 in range(<np.uint64_t>self.s1):
                mm_s0._set_coarse(mi1)
            mm_s0._compress(mm_s)
            return
        else:
            self.recursive_morton_mask(level, pos, dds, mi1, cur_ind)
        # Set coarse morton indices in order
        self.set_coarse_bool(mm_s0, mm_g0)
        self.set_refined_list(mm_s0, mm_g0)
        self.set_refined_bool(mm_s0, mm_g0)
        # Compress
        mm_s0._compress(mm_s)
        mm_g0._compress(mm_g)

    def find_files(self,
                   np.ndarray[np.uint8_t, ndim=1] file_mask_p,
                   np.ndarray[np.uint8_t, ndim=1] file_mask_g):
        cdef np.uint64_t i
        cdef np.int32_t level = 0
        cdef np.uint64_t mi1
        mi1 = ~(<np.uint64_t>0)
        cdef np.float64_t pos[3]
        cdef np.float64_t dds[3]
        for i in range(3):
            pos[i] = self.DLE[i]
            dds[i] = self.DRE[i] - self.DLE[i]
        # Fill with input
        for i in range(self.nfiles):
            self.file_mask_p[i] = file_mask_p[i]
            self.file_mask_g[i] = file_mask_g[i]
        # Recurse
        self.recursive_morton_files(level, pos, dds, mi1)
        # Fill with results
        for i in range(self.nfiles):
            file_mask_p[i] = self.file_mask_p[i]
            if file_mask_p[i]:
                file_mask_g[i] = 0
            else:
                file_mask_g[i] = self.file_mask_g[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef bint is_refined(self, np.uint64_t mi1):
        return self.bitmap.collisions._isref(mi1)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef bint is_refined_files(self, np.uint64_t mi1):
        cdef np.uint64_t i
        if self.bitmap.collisions._isref(mi1):
            # Don't refine if files all selected already
            for i in range(self.nfiles):
                if self.file_mask_p[i] == 0:
                    if self.bitmap.bitmasks._isref(i, mi1) == 1:
                        return 1
            return 0
        else:
            return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void add_coarse(self, np.uint64_t mi1, int bbox = 2):
        self.coarse_select_bool[mi1] = 1
        # Neighbors
        if (self.ngz > 0) and (bbox == 2):
            if not self.is_refined(mi1):
                self.add_neighbors_coarse(mi1)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void set_files_coarse(self, np.uint64_t mi1):
        cdef np.uint64_t i
        cdef bint flag_ref = self.is_refined(mi1)
        # Flag files at coarse level
        if flag_ref == 0:
            for i in range(self.nfiles):
                if self.file_mask_p[i] == 0:
                    if self.bitmap.bitmasks._get_coarse(i, mi1) == 1:
                        self.file_mask_p[i] = 1
        # Neighbors
        if (flag_ref == 0) and (self.ngz > 0):
            self.set_files_neighbors_coarse(mi1)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int add_refined(self, np.uint64_t mi1, np.uint64_t mi2, int bbox = 2) except -1:
        self.refined_select_bool[mi2] = 1
        # Neighbors
        if (self.ngz > 0) and (bbox == 2):
            self.add_neighbors_refined(mi1, mi2)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void set_files_refined(self, np.uint64_t mi1, np.uint64_t mi2):
        cdef np.uint64_t i
        # Flag files
        for i in range(self.nfiles):
            if self.file_mask_p[i] == 0:
                if self.bitmap.bitmasks._get(i, mi1, mi2):
                    self.file_mask_p[i] = 1
        # Neighbors
        if (self.ngz > 0):
            self.set_files_neighbors_refined(mi1, mi2)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void add_neighbors_coarse(self, np.uint64_t mi1):
        cdef np.uint64_t m
        cdef np.uint32_t ntot
        cdef np.uint64_t mi1_n
        ntot = morton_neighbors_coarse(mi1, self.max_index1,
                                       self.periodicity,
                                       self.ngz, self.neighbors,
                                       self.ind1_n, self.neighbor_list1)
        for m in range(ntot):
            mi1_n = self.neighbor_list1[m]
            self.coarse_ghosts_bool[mi1_n] = 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void set_files_neighbors_coarse(self, np.uint64_t mi1):
        cdef np.uint64_t i, m
        cdef np.uint32_t ntot
        cdef np.uint64_t mi1_n
        ntot = morton_neighbors_coarse(mi1, self.max_index1,
                                       self.periodicity,
                                       self.ngz, self.neighbors,
                                       self.ind1_n, self.neighbor_list1)
        for m in range(ntot):
            mi1_n = self.neighbor_list1[m]
            for i in range(self.nfiles):
                if self.file_mask_g[i] == 0:
                    if self.bitmap.bitmasks._get_coarse(i, mi1_n):
                        self.file_mask_g[i] = 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void add_neighbors_refined(self, np.uint64_t mi1, np.uint64_t mi2):
        cdef int m
        cdef np.uint32_t ntot
        cdef np.uint64_t mi1_n, mi2_n
        ntot = morton_neighbors_refined(mi1, mi2,
                                        self.max_index1, self.max_index2,
                                        self.periodicity, self.ngz,
                                        self.neighbors, self.ind1_n, self.ind2_n,
                                        self.neighbor_list1, self.neighbor_list2)
        for m in range(<np.int32_t>ntot):
            mi1_n = self.neighbor_list1[m]
            mi2_n = self.neighbor_list2[m]
            self.coarse_ghosts_bool[mi1_n] = 1
            IF RefinedExternalGhosts == 1:
                if mi1_n == mi1:
                    self.refined_ghosts_bool[mi2_n] = 1
                else:
                    self.refined_ghosts_list._set(mi1_n, mi2_n)
            ELSE:
                if mi1_n == mi1:
                    self.refined_ghosts_bool[mi2_n] = 1
                elif self.is_refined(mi1_n) == 1:
                    self.refined_ghosts_list._set(mi1_n, mi2_n)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void set_files_neighbors_refined(self, np.uint64_t mi1, np.uint64_t mi2):
        cdef int i, m
        cdef np.uint32_t ntot
        cdef np.uint64_t mi1_n, mi2_n
        ntot = morton_neighbors_refined(mi1, mi2,
                                        self.max_index1, self.max_index2,
                                        self.periodicity, self.ngz,
                                        self.neighbors, self.ind1_n, self.ind2_n,
                                        self.neighbor_list1, self.neighbor_list2)
        for m in range(<np.int32_t>ntot):
            mi1_n = self.neighbor_list1[m]
            mi2_n = self.neighbor_list2[m]
            if self.is_refined(mi1_n) == 1:
                for i in range(self.nfiles):
                    if self.file_mask_g[i] == 0:
                        if self.bitmap.bitmasks._get(i, mi1_n, mi2_n) == 1:
                            self.file_mask_g[i] = 1
            else:
                for i in range(self.nfiles):
                    if self.file_mask_g[i] == 0:
                        if self.bitmap.bitmasks._get_coarse(i, mi1_n) == 1:
                            self.file_mask_g[i] = 1
                            break # If not refined, only one file should be selected

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void set_coarse_list(self, BoolArrayColl mm_s, BoolArrayColl mm_g):
        self.coarse_select_list._fill_bool(mm_s)
        self.coarse_ghosts_list._fill_bool(mm_g)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void set_refined_list(self, BoolArrayColl mm_s, BoolArrayColl mm_g):
        self.refined_ghosts_list._fill_bool(mm_g)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void set_coarse_bool(self, BoolArrayColl mm_s, BoolArrayColl mm_g):
        cdef np.uint64_t mi1
        mm_s._set_coarse_array_ptr(&self.coarse_select_bool[0])
        for mi1 in range(self.s1):
            self.coarse_select_bool[mi1] = 0
        mm_g._set_coarse_array_ptr(&self.coarse_ghosts_bool[0])
        for mi1 in range(self.s1):
            self.coarse_ghosts_bool[mi1] = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void set_refined_bool(self, BoolArrayColl mm_s, BoolArrayColl mm_g):
        mm_s._append(self.select_ewah)
        mm_g._append(self.ghosts_ewah)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void push_refined_bool(self, np.uint64_t mi1):
        cdef np.uint64_t mi2
        self.select_ewah._set_refined_array_ptr(mi1, &self.refined_select_bool[0])
        for mi2 in range(self.s2):
            self.refined_select_bool[mi2] = 0
        self.ghosts_ewah._set_refined_array_ptr(mi1, &self.refined_ghosts_bool[0])
        for mi2 in range(self.s2):
            self.refined_ghosts_bool[mi2] = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void add_ghost_zones(self, BoolArrayColl mm_s, BoolArrayColl mm_g):
        cdef np.uint64_t mi1, mi2, mi1_n, mi2_n
        # Get ghost zones, unordered
        for mi1 in range(self.s1):
            if mm_s._get_coarse(mi1):
                if self.is_refined(mi1):
                    for mi2 in range(self.s2):
                        if mm_s._get(mi1, mi2):
                            self.add_neighbors_refined(mi1, mi2)
                    # self.push_refined_bool(mi1)
                    self.ghosts_ewah._set_refined_array_ptr(mi1,
                                                            &self.refined_ghosts_bool[0])
                    for mi2 in range(self.s2):
                        self.refined_ghosts_bool[mi2] = 0
                else:
                    self.add_neighbors_coarse(mi1)
        # Add ghost zones to bool array in order
        mm_g._set_coarse_array_ptr(&self.coarse_ghosts_bool[0])
        for mi1 in range(self.s1):
            self.coarse_ghosts_bool[mi1] = 0
        self.refined_ghosts_list._fill_bool(mm_g)
        mm_g._append(self.ghosts_ewah)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int fill_subcells_mi1(self,
                               np.uint64_t nlevel,
                               np.uint64_t ind1[3]) except -1:
        cdef np.uint64_t imi, fmi
        cdef np.uint64_t mi
        cdef np.uint64_t shift_by = 3 * (self.bitmap.index_order1 - nlevel)
        imi = encode_morton_64bit(ind1[0], ind1[1], ind1[2]) << shift_by
        fmi = imi + (1 << shift_by)
        for mi in range(imi, fmi):
            self.add_coarse(mi, 1)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int fill_subcells_mi2(self,
                               np.uint64_t nlevel,
                               np.uint64_t mi1,
                               np.uint64_t ind2[3]) except -1:
        cdef np.uint64_t imi, fmi
        cdef np.uint64_t shift_by = 3 * ((self.bitmap.index_order2 +
                                          self.bitmap.index_order1) - nlevel)
        imi = encode_morton_64bit(ind2[0], ind2[1], ind2[2]) << shift_by
        fmi = imi + (1 << shift_by)
        for mi2 in range(imi, fmi):
            self.add_refined(mi1, mi2, 1)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int recursive_morton_mask(
        self, np.int32_t level, np.float64_t pos[3],
        np.float64_t dds[3], np.uint64_t mi1, np.uint64_t cur_ind[3]) except -1:
        cdef np.uint64_t mi2
        cdef np.float64_t npos[3]
        cdef np.float64_t rpos[3]
        cdef np.float64_t ndds[3]
        cdef np.uint64_t nlevel
        cdef np.uint64_t ind1[3]
        cdef np.uint64_t ind2[3]
        cdef np.uint64_t ncur_ind[3]
        cdef np.uint64_t* zeros = [0, 0, 0]
        cdef int i, j, k, m, sbbox
        PyErr_CheckSignals()
        for i in range(3):
            ndds[i] = dds[i]/2
        nlevel = level + 1
        # Loop over octs
        for i in range(2):
            npos[0] = pos[0] + i*ndds[0]
            rpos[0] = npos[0] + ndds[0]
            ncur_ind[0] = (cur_ind[0] << 1) + i
            for j in range(2):
                npos[1] = pos[1] + j*ndds[1]
                rpos[1] = npos[1] + ndds[1]
                ncur_ind[1] = (cur_ind[1] << 1) + j
                for k in range(2):
                    npos[2] = pos[2] + k*ndds[2]
                    rpos[2] = npos[2] + ndds[2]
                    ncur_ind[2] = (cur_ind[2] << 1) + k
                    # Only recurse into selected cells
                    sbbox = self.selector.select_bbox_edge(npos, rpos)
                    if sbbox == 0:
                        continue
                    if nlevel < self.order1:
                        if sbbox == 1:
                            self.fill_subcells_mi1(nlevel, ncur_ind)
                        else:
                            self.recursive_morton_mask(
                                nlevel, npos, ndds, mi1, ncur_ind)
                    elif nlevel == self.order1:
                        mi1 = encode_morton_64bit(
                            ncur_ind[0], ncur_ind[1], ncur_ind[2])
                        if sbbox == 2: # an edge cell
                            if self.is_refined(mi1) == 1:
                                # note we pass zeros here in the last argument
                                # this is because we now need to generate
                                # *refined* indices above order1 so we need to
                                # start a new running count of refined indices.
                                #
                                # note that recursive_morton_mask does not
                                # mutate the last argument (a new index is
                                # calculated in each stack frame) so this is
                                # safe
                                self.recursive_morton_mask(
                                    nlevel, npos, ndds, mi1, zeros)
                        self.add_coarse(mi1, sbbox)
                        self.push_refined_bool(mi1)
                    elif nlevel < (self.order1 + self.order2):
                        if sbbox == 1:
                            self.fill_subcells_mi2(nlevel, mi1, ncur_ind)
                        else:
                            self.recursive_morton_mask(
                                nlevel, npos, ndds, mi1, ncur_ind)
                    elif nlevel == (self.order1 + self.order2):
                        mi2 = encode_morton_64bit(
                            ncur_ind[0], ncur_ind[1], ncur_ind[2])
                        self.add_refined(mi1, mi2, sbbox)
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void recursive_morton_files(self, np.int32_t level, np.float64_t pos[3],
                                     np.float64_t dds[3], np.uint64_t mi1):
        cdef np.uint64_t mi2
        cdef np.float64_t npos[3]
        cdef np.float64_t rpos[3]
        cdef np.float64_t ndds[3]
        cdef np.uint64_t nlevel
        cdef np.float64_t DLE[3]
        cdef np.uint64_t ind1[3]
        cdef np.uint64_t ind2[3]
        cdef int i, j, k, m
        for i in range(3):
            ndds[i] = dds[i]/2
        nlevel = level + 1
        # Loop over octs
        for i in range(2):
            npos[0] = pos[0] + i*ndds[0]
            rpos[0] = npos[0] + ndds[0]
            for j in range(2):
                npos[1] = pos[1] + j*ndds[1]
                rpos[1] = npos[1] + ndds[1]
                for k in range(2):
                    npos[2] = pos[2] + k*ndds[2]
                    rpos[2] = npos[2] + ndds[2]
                    # Only recurse into selected cells
                    if not self.selector.select_bbox(npos, rpos): continue
                    if nlevel < self.order1:
                        self.recursive_morton_files(nlevel, npos, ndds, mi1)
                    elif nlevel == self.order1:
                        mi1 = bounded_morton_dds(npos[0], npos[1], npos[2], self.DLE, ndds)
                        if self.is_refined_files(mi1):
                            self.recursive_morton_files(nlevel, npos, ndds, mi1)
                        self.set_files_coarse(mi1)
                    elif nlevel < (self.order1 + self.order2):
                        self.recursive_morton_files(nlevel, npos, ndds, mi1)
                    elif nlevel == (self.order1 + self.order2):
                        decode_morton_64bit(mi1,ind1)
                        for m in range(3):
                            DLE[m] = self.DLE[m] + ndds[m]*ind1[m]*self.max_index2
                        mi2 = bounded_morton_dds(npos[0], npos[1], npos[2], DLE, ndds)
                        self.set_files_refined(mi1,mi2)

cdef class ParticleBitmapOctreeContainer(SparseOctreeContainer):
    cdef Oct** oct_list
    cdef public int max_level
    cdef public int n_ref
    cdef int loaded # Loaded with load_octree?
    cdef np.uint8_t* _ptr_index_base_roots
    cdef np.uint8_t* _ptr_index_base_octs
    cdef np.uint64_t* _ptr_octs_per_root
    cdef public np.uint8_t[:] _index_base_roots
    cdef public np.uint8_t[:] _index_base_octs
    cdef np.uint64_t[:] _octs_per_root
    cdef public int overlap_cells
    def __init__(self, domain_dimensions, domain_left_edge, domain_right_edge,
                 int num_root, over_refine = 1):
        super(ParticleBitmapOctreeContainer, self).__init__(
            domain_dimensions, domain_left_edge, domain_right_edge,
            over_refine)
        self.loaded = 0
        self.fill_style = "o"
        self.partial_coverage = 2
        self.overlap_cells = 0
        # Now the overrides
        self.max_level = -1
        self.max_root = num_root
        self.root_nodes = <OctKey*> malloc(sizeof(OctKey) * num_root)
        self._ptr_index_base_roots = <np.uint8_t*> malloc(sizeof(np.uint8_t) * num_root)
        self._ptr_octs_per_root = <np.uint64_t*> malloc(sizeof(np.uint64_t) * num_root)
        for i in range(num_root):
            self.root_nodes[i].key = -1
            self.root_nodes[i].node = NULL
            self._ptr_index_base_roots[i] = 1
            self._ptr_octs_per_root[i] = 0
        self._index_base_roots = <np.uint8_t[:num_root]> self._ptr_index_base_roots
        self._octs_per_root = <np.uint64_t[:num_root]> self._ptr_octs_per_root

    def allocate_domains(self, counts = None):
        if counts is None:
            counts = [self.max_root]
        OctreeContainer.allocate_domains(self, counts)

    def finalize(self):
        # Assign domain ind
        cdef SelectorObject selector = AlwaysSelector(None)
        selector.overlap_cells = self.overlap_cells
        cdef oct_visitors.AssignDomainInd visitor
        visitor = oct_visitors.AssignDomainInd(self)
        self.visit_all_octs(selector, visitor)
        assert ((visitor.global_index+1)*visitor.nz == visitor.index)
        # Copy indexes
        self._ptr_index_base_octs = <np.uint8_t*> malloc(sizeof(np.uint8_t)*self.nocts)
        self._index_base_octs = <np.uint8_t[:self.nocts]> self._ptr_index_base_octs
        cdef np.int64_t nprev_octs = 0
        cdef int i
        for i in range(self.num_root):
            self._index_base_octs[nprev_octs:(nprev_octs+self._octs_per_root[i])] = self._index_base_roots[i]
            nprev_octs += self._octs_per_root[i]

    cdef visit_assign(self, Oct *o, np.int64_t *lpos, int level, int *max_level,
                      np.int64_t index_root):
        cdef int i, j, k
        if o.children == NULL:
            self.oct_list[lpos[0]] = o
            self._index_base_octs[lpos[0]] = self._index_base_roots[index_root]
            lpos[0] += 1
        max_level[0] = imax(max_level[0], level)
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if o.children != NULL \
                       and o.children[cind(i,j,k)] != NULL:
                        self.visit_assign(o.children[cind(i,j,k)], lpos,
                                          level + 1, max_level, index_root)
        return

    cdef Oct* allocate_oct(self):
        #Allocate the memory, set to NULL or -1
        #We reserve space for n_ref particles, but keep
        #track of how many are used with np initially 0
        self.nocts += 1
        cdef Oct *my_oct = <Oct*> malloc(sizeof(Oct))
        my_oct.domain = -1
        my_oct.file_ind = 0
        my_oct.domain_ind = self.nocts - 1
        my_oct.children = NULL
        return my_oct

    def get_index_base_octs(self, np.int64_t[:] domain_ind):
        cdef np.int64_t ndst = np.max(domain_ind) + 1
        ind = np.zeros(ndst, 'int64') - 1
        self._get_index_base_octs(ind, domain_ind)
        return ind[ind >= 0]

    cdef void _get_index_base_octs(self, np.int64_t[:] ind, np.int64_t[:] domain_ind):
        cdef SelectorObject selector = AlwaysSelector(None)
        selector.overlap_cells = self.overlap_cells
        cdef oct_visitors.IndexMaskMapOcts visitor
        visitor = oct_visitors.IndexMaskMapOcts(self)
        visitor.oct_mask = self._index_base_octs
        visitor.oct_index = ind
        visitor.map_domain_ind = domain_ind
        self.visit_all_octs(selector, visitor)

    def __dealloc__(self):
        #Call the freemem ops on every ocy
        #of the root mesh recursively
        cdef int i
        if self.root_nodes== NULL: return
        if self.loaded == 0:
            for i in range(self.max_root):
                if self.root_nodes[i].node == NULL: continue
                self.visit_free(&self.root_nodes.node[i], 0)
            self.root_nodes = NULL
        free(self.oct_list)
        free(self._ptr_index_base_roots)
        free(self._ptr_index_base_octs)
        free(self._ptr_octs_per_root)
        self.oct_list = NULL

    cdef void visit_free(self, Oct *o, int free_this):
        #Free the memory for this oct recursively
        cdef int i, j, k
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if o.children != NULL \
                       and o.children[cind(i,j,k)] != NULL:
                        self.visit_free(o.children[cind(i,j,k)], 1)
        if o.children != NULL:
            free(o.children)
        if free_this == 1:
            free(o)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void recursive_add(self, Oct *o, np.ndarray[np.uint64_t, ndim=1] indices,
                            int level, int *max_level, int domain_id, int *count):
        cdef np.int64_t no = indices.shape[0], beg, end, nind
        cdef np.int64_t index
        cdef int i, j, k
        cdef int ind[3]
        cdef Oct *noct
        cdef Oct *noct_ch
        beg = end = 0
        if level > max_level[0]: max_level[0] = level
        # Initialize children
        if o.children == NULL:
            o.children = <Oct **> malloc(sizeof(Oct *)*8)
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        o.children[cind(i,j,k)] = NULL
                        # noct = self.allocate_oct()
                        # noct.domain = o.domain
                        # noct.file_ind = 0
                        # o.children[cind(i,j,k)] = noct
        # Loop through sets of particles with matching prefix at this level
        while end < no:
            beg = end
            index = (indices[beg] >> ((ORDER_MAX - level)*3))
            while (end < no) and (index == (indices[end] >> ((ORDER_MAX - level)*3))):
                end += 1
            nind = (end - beg)
            # Add oct
            for i in range(3):
                ind[i] = ((index >> (2 - i)) & 1)
            # noct = o.children[cind(ind[0],ind[1],ind[2])]
            if o.children[cind(ind[0],ind[1],ind[2])] != NULL:
                raise Exception('Child was already initialized...')
            noct = self.allocate_oct()
            noct.domain = o.domain
            o.children[cind(ind[0],ind[1],ind[2])] = noct
            # Don't add it to the list if it will be refined
            if nind > self.n_ref and level < ORDER_MAX:
                self.nocts -= 1
                noct.domain_ind = -1 # overwritten by finalize
            else:
                count[0] += 1
            noct.file_ind = o.file_ind
            # noct.file_ind = nind
            # o.file_ind = self.n_ref + 1
            # Refine oct or add its children
            if nind > self.n_ref and level < ORDER_MAX:
                self.recursive_add(noct, indices[beg:end], level+1,
                                   max_level, domain_id, count)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def add(self, np.ndarray[np.uint64_t, ndim=1] indices,
             np.uint64_t order1, int domain_id = -1):
        #Add this particle to the root oct
        #Then if that oct has children, add it to them recursively
        #If the child needs to be refined because of max particles, do so
        cdef Oct *root = NULL
        cdef np.int64_t no = indices.shape[0], beg, end, index, nind
        cdef int i, level
        cdef int ind[3]
        cdef np.uint64_t ind64[3]
        cdef int max_level = self.max_level
        # Note what we're doing here: we have decided the root will always be
        # zero, since we're in a forest of octrees, where the root_mesh node is
        # the level 0.  This means our morton indices should be made with
        # respect to that, which means we need to keep a few different arrays
        # of them.
        cdef np.int64_t index_root = 0
        cdef int root_count
        beg = end = 0
        self._octs_per_root[:] = 1 # Roots count reguardless
        while end < no:
            # Determine number of octs with this prefix
            beg = end
            index = (indices[beg] >> ((ORDER_MAX - self.level_offset)*3))
            while (end < no) and (index == (indices[end] >> ((ORDER_MAX - self.level_offset)*3))):
                end += 1
            nind = (end - beg)
            # Find root for prefix
            decode_morton_64bit(index, ind64)
            for i in range(3):
                ind[i] = ind64[i]
            while (index_root < self.num_root) and \
                  (self.ipos_to_key(ind) != self.root_nodes[index_root].key):
                index_root += 1
            if index_root >= self.num_root:
                raise Exception('No root found for {},{},{}'.format(ind[0],ind[1],ind[2]))
            root = self.root_nodes[index_root].node
            # self.get_root(ind, &root)
            # if root == NULL:
            #     raise Exception('No root found for {},{},{}'.format(ind[0],ind[1],ind[2]))
            root.file_ind = index_root
            # Refine root as necessary
            if (end - beg) > self.n_ref:
                root_count = 0
                self.nocts -= 1
                self.recursive_add(root, indices[beg:end], self.level_offset+1,
                                   &max_level, domain_id, &root_count)
                self._octs_per_root[index_root] = <np.uint64_t>root_count
        self.max_level = max_level
        assert(self.nocts == np.sum(self._octs_per_root))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Oct *refine_oct(self, Oct *o, np.uint64_t index, int level):
        #Allocate and initialize child octs
        #Attach particles to child octs
        #Remove particles from this oct entirely
        cdef int i, j, k
        cdef int ind[3]
        cdef Oct *noct

        # Initialize empty children
        if o.children == NULL:
            o.children = <Oct **> malloc(sizeof(Oct *)*8)

        # This version can be used to just add the child containing the index
        #     for i in range(2):
        #         for j in range(2):
        #             for k in range(2):
        #                 o.children[cind(i,j,k)] = NULL
        # # Only allocate and count the indexed oct
        # for i in range(3):
        #     ind[i] = (index >> ((ORDER_MAX - level)*3 + (2 - i))) & 1

        # noct = self.allocate_oct()
        # noct.domain = o.domain
        # noct.file_ind = 0
        # o.children[cind(ind[0],ind[1],ind[2])] = noct
        # o.file_ind = self.n_ref + 1


        for i in range(2):
            for j in range(2):
                for k in range(2):
                    noct = self.allocate_oct()
                    noct.domain = o.domain
                    noct.file_ind = 0
                    o.children[cind(i,j,k)] = noct
        o.file_ind = self.n_ref + 1
        for i in range(3):
            ind[i] = (index >> ((ORDER_MAX - level)*3 + (2 - i))) & 1
        noct = o.children[cind(ind[0],ind[1],ind[2])]
        return noct

    cdef void filter_particles(self, Oct *o, np.uint64_t *data, np.int64_t p,
                               int level):
        # Now we look at the last nref particles to decide where they go.
        cdef int n = imin(p, self.n_ref)
        cdef np.uint64_t *arr = data + imax(p - self.n_ref, 0)
        # Now we figure out our prefix, which is the oct address at this level.
        # As long as we're actually in Morton order, we do not need to worry
        # about *any* of the other children of the oct.
        prefix1 = data[p] >> (ORDER_MAX - level)*3
        for i in range(n):
            prefix2 = arr[i] >> (ORDER_MAX - level)*3
            if (prefix1 == prefix2):
                o.file_ind += 1
