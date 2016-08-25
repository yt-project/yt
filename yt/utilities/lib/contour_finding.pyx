"""
A two-pass contour finding algorithm



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free, realloc
from yt.geometry.selection_routines cimport \
    SelectorObject, AlwaysSelector, OctreeSubsetSelector, \
    anyfloat
from yt.utilities.lib.fp_utils cimport imax
from yt.geometry.oct_container cimport \
    OctreeContainer, OctInfo
from yt.geometry.oct_visitors cimport \
    Oct
from .amr_kdtools cimport Node
from .partitioned_grid cimport \
    PartitionedGrid
from .volume_container cimport \
    VolumeContainer, vc_index, vc_pos_index
import sys

cdef inline ContourID *contour_create(np.int64_t contour_id,
                               ContourID *prev = NULL):
    node = <ContourID *> malloc(sizeof(ContourID))
    #print "Creating contour with id", contour_id
    node.contour_id = contour_id
    node.next = node.parent = NULL
    node.prev = prev
    node.count = 1
    if prev != NULL: prev.next = node
    return node

cdef inline void contour_delete(ContourID *node):
    if node.prev != NULL: node.prev.next = node.next
    if node.next != NULL: node.next.prev = node.prev
    free(node)

cdef inline ContourID *contour_find(ContourID *node):
    cdef ContourID *temp
    cdef ContourID *root
    root = node
    # First we find the root
    while root.parent != NULL and root.parent != root:
        root = root.parent
    if root == root.parent: root.parent = NULL
    # Now, we update everything along the tree.
    # So now everything along the line to the root has the parent set to the
    # root.
    while node.parent != NULL:
        temp = node.parent
        root.count += node.count
        node.count = 0
        node.parent = root
        node = temp
    return root

cdef inline void contour_union(ContourID *node1, ContourID *node2):
    if node1 == node2:
        return
    node1 = contour_find(node1)
    node2 = contour_find(node2)
    if node1 == node2:
        return
    cdef ContourID *pri
    cdef ContourID *sec
    if node1.count > node2.count:
        pri = node1
        sec = node2
    elif node2.count > node1.count:
        pri = node2
        sec = node1
    # might be a tie
    elif node1.contour_id < node2.contour_id:
        pri = node1
        sec = node2
    else:
        pri = node2
        sec = node1
    pri.count += sec.count
    sec.count = 0
    sec.parent = pri

cdef inline int candidate_contains(CandidateContour *first,
                            np.int64_t contour_id,
                            np.int64_t join_id = -1):
    while first != NULL:
        if first.contour_id == contour_id \
            and first.join_id == join_id: return 1
        first = first.next
    return 0

cdef inline CandidateContour *candidate_add(CandidateContour *first,
                                     np.int64_t contour_id,
                                     np.int64_t join_id = -1):
    cdef CandidateContour *node
    node = <CandidateContour *> malloc(sizeof(CandidateContour))
    node.contour_id = contour_id
    node.join_id = join_id
    node.next = first
    return node

cdef class ContourTree:
    # This class is essentially a Union-Find algorithm.  What we want to do is
    # to, given a connection between two objects, identify the unique ID for
    # those two objects.  So what we have is a collection of contours, and they
    # eventually all get joined and contain lots of individual IDs.  But it's
    # easy to find the *first* contour, i.e., the primary ID, for each of the
    # subsequent IDs.
    #
    # This means that we can connect id 202483 to id 2472, and if id 2472 is
    # connected to id 143, the connection will *actually* be from 202483 to
    # 143.  In this way we can speed up joining things and knowing their
    # "canonical" id.
    #
    # This is a multi-step process, since we first want to connect all of the
    # contours, then we end up wanting to coalesce them, and ultimately we join
    # them at the end.  The join produces a table that maps the initial to the
    # final, and we can go through and just update all of those.
    cdef ContourID *first
    cdef ContourID *last

    def clear(self):
        # Here, we wipe out ALL of our contours, but not the pointers to them
        cdef ContourID *cur
        cdef ContourID *next
        cur = self.first
        while cur != NULL:
            next = cur.next
            free(cur)
            cur = next
        self.first = self.last = NULL

    def __init__(self):
        self.first = self.last = NULL

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def add_contours(self, np.ndarray[np.int64_t, ndim=1] contour_ids):
        # This adds new contours, from the given contour IDs, to the tree.
        # Each one can be connected to a parent, as well as to next/prev in the
        # set of contours belonging to this tree.
        cdef int i, n
        n = contour_ids.shape[0]
        cdef ContourID *cur = self.last
        for i in range(n):
            #print i, contour_ids[i]
            cur = contour_create(contour_ids[i], cur)
            if self.first == NULL: self.first = cur
        self.last = cur

    def add_contour(self, np.int64_t contour_id):
        self.last = contour_create(contour_id, self.last)

    def cull_candidates(self, np.ndarray[np.int64_t, ndim=3] candidates):
        # This function looks at each preliminary contour ID belonging to a
        # given collection of values, and then if need be it creates a new
        # contour for it.
        cdef int i, j, k, ni, nj, nk, nc
        cdef CandidateContour *first = NULL
        cdef CandidateContour *temp
        cdef np.int64_t cid
        nc = 0
        ni = candidates.shape[0]
        nj = candidates.shape[1]
        nk = candidates.shape[2]
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    cid = candidates[i,j,k]
                    if cid == -1: continue
                    if candidate_contains(first, cid) == 0:
                        nc += 1
                        first = candidate_add(first, cid)
        cdef np.ndarray[np.int64_t, ndim=1] contours
        contours = np.empty(nc, dtype="int64")
        i = 0
        # This removes all the temporary contours for this set of contours and
        # instead constructs a final list of them.
        while first != NULL:
            contours[i] = first.contour_id
            i += 1
            temp = first.next
            free(first)
            first = temp
        return contours

    def cull_joins(self, np.ndarray[np.int64_t, ndim=2] cjoins):
        # This coalesces contour IDs, so that we have only the final name
        # resolutions -- the .join_id from a candidate.  So many items will map
        # to a single join_id.
        cdef int i, ni, nc
        cdef CandidateContour *first = NULL
        cdef CandidateContour *temp
        cdef np.int64_t cid1, cid2
        nc = 0
        ni = cjoins.shape[0]
        for i in range(ni):
            cid1 = cjoins[i,0]
            cid2 = cjoins[i,1]
            if cid1 == -1: continue
            if cid2 == -1: continue
            if candidate_contains(first, cid1, cid2) == 0:
                nc += 1
                first = candidate_add(first, cid1, cid2)
        cdef np.ndarray[np.int64_t, ndim=2] contours
        contours = np.empty((nc,2), dtype="int64")
        i = 0
        while first != NULL:
            contours[i,0] = first.contour_id
            contours[i,1] = first.join_id
            i += 1
            temp = first.next
            free(first)
            first = temp
        return contours

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def add_joins(self, np.ndarray[np.int64_t, ndim=2] join_tree):
        cdef int i, n, ins
        cdef np.int64_t cid1, cid2
        # Okay, this requires lots of iteration, unfortunately
        cdef ContourID *cur
        cdef ContourID *c1
        cdef ContourID *c2
        n = join_tree.shape[0]
        #print "Counting"
        #print "Checking", self.count()
        for i in range(n):
            ins = 0
            cid1 = join_tree[i, 0]
            cid2 = join_tree[i, 1]
            c1 = c2 = NULL
            cur = self.first
            #print "Looking for ", cid1, cid2
            while c1 == NULL or c2 == NULL:
                if cur.contour_id == cid1:
                    c1 = contour_find(cur)
                if cur.contour_id == cid2:
                    c2 = contour_find(cur)
                ins += 1
                cur = cur.next
                if cur == NULL: break
            if c1 == NULL or c2 == NULL:
                if c1 == NULL: print "  Couldn't find ", cid1
                if c2 == NULL: print "  Couldn't find ", cid2
                print "  Inspected ", ins
                raise RuntimeError
            else:
                c1.count = c2.count = 0
                contour_union(c1, c2)

    def count(self):
        cdef int n = 0
        cdef ContourID *cur = self.first
        while cur != NULL:
            cur = cur.next
            n += 1
        return n

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def export(self):
        cdef int n = self.count()
        cdef ContourID *cur
        cdef ContourID *root
        cur = self.first
        cdef np.ndarray[np.int64_t, ndim=2] joins
        joins = np.empty((n, 2), dtype="int64")
        n = 0
        while cur != NULL:
            root = contour_find(cur)
            joins[n, 0] = cur.contour_id
            joins[n, 1] = root.contour_id
            cur = cur.next
            n += 1
        return joins

    def __dealloc__(self):
        self.clear()

cdef class TileContourTree:
    cdef np.float64_t min_val
    cdef np.float64_t max_val

    def __init__(self, np.float64_t min_val, np.float64_t max_val):
        self.min_val = min_val
        self.max_val = max_val

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def identify_contours(self, np.ndarray[np.float64_t, ndim=3] values,
                                np.ndarray[np.int64_t, ndim=3] contour_ids,
                                np.ndarray[np.uint8_t, ndim=3] mask,
                                np.int64_t start):
        # This just looks at neighbor values and tries to identify which zones
        # are touching by face within a given brick.
        cdef int i, j, k, ni, nj, nk, offset
        cdef int off_i, off_j, off_k, oi, ok, oj
        cdef ContourID *cur = NULL
        cdef ContourID *c1
        cdef ContourID *c2
        cdef np.float64_t v
        cdef np.int64_t nc
        ni = values.shape[0]
        nj = values.shape[1]
        nk = values.shape[2]
        nc = 0
        cdef ContourID **container = <ContourID**> malloc(
                sizeof(ContourID*)*ni*nj*nk)
        for i in range(ni*nj*nk): container[i] = NULL
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    v = values[i,j,k]
                    if mask[i,j,k] == 0: continue
                    if v < self.min_val or v > self.max_val: continue
                    nc += 1
                    c1 = contour_create(nc + start)
                    cur = container[i*nj*nk + j*nk + k] = c1
                    for oi in range(3):
                        off_i = oi - 1 + i
                        if not (0 <= off_i < ni): continue
                        for oj in range(3):
                            off_j = oj - 1 + j
                            if not (0 <= off_j < nj): continue
                            for ok in range(3):
                                if oi == oj == ok == 1: continue
                                off_k = ok - 1 + k
                                if not (0 <= off_k < nk): continue
                                if off_k > k and off_j > j and off_i > i:
                                    continue
                                offset = off_i*nj*nk + off_j*nk + off_k
                                c2 = container[offset]
                                if c2 == NULL: continue
                                c2 = contour_find(c2)
                                cur.count = c2.count = 0
                                contour_union(cur, c2)
                                cur = contour_find(cur)
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    c1 = container[i*nj*nk + j*nk + k]
                    if c1 == NULL: continue
                    c1 = contour_find(c1)
                    contour_ids[i,j,k] = c1.contour_id

        for i in range(ni*nj*nk):
            if container[i] != NULL: free(container[i])
        free(container)
        return nc

@cython.boundscheck(False)
@cython.wraparound(False)
def link_node_contours(Node trunk, contours, ContourTree tree,
        np.ndarray[np.int64_t, ndim=1] node_ids):
    cdef int n_nodes = node_ids.shape[0]
    cdef np.int64_t node_ind
    cdef VolumeContainer **vcs = <VolumeContainer **> malloc(
        sizeof(VolumeContainer*) * n_nodes)
    cdef int i
    cdef PartitionedGrid pg
    for i in range(n_nodes):
        pg = contours[node_ids[i]][2]
        vcs[i] = pg.container
    cdef np.ndarray[np.uint8_t] examined = np.zeros(n_nodes, "uint8")
    for _, cinfo in sorted(contours.items(), key = lambda a: -a[1][0]):
        _, node_ind, pg, _ = cinfo
        construct_boundary_relationships(trunk, tree, node_ind,
            examined, vcs, node_ids)
        examined[node_ind] = 1

cdef inline void get_spos(VolumeContainer *vc, int i, int j, int k,
                          int axis, np.float64_t *spos):
    spos[0] = vc.left_edge[0] + i * vc.dds[0]
    spos[1] = vc.left_edge[1] + j * vc.dds[1]
    spos[2] = vc.left_edge[2] + k * vc.dds[2]
    spos[axis] += 0.5 * vc.dds[axis]

cdef inline int spos_contained(VolumeContainer *vc, np.float64_t *spos):
    cdef int i
    for i in range(3):
        if spos[i] <= vc.left_edge[i] or spos[i] >= vc.right_edge[i]: return 0
    return 1

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void construct_boundary_relationships(Node trunk, ContourTree tree,
                np.int64_t nid, np.ndarray[np.uint8_t, ndim=1] examined,
                VolumeContainer **vcs,
                np.ndarray[np.int64_t, ndim=1] node_ids):
    # We only look at the boundary and find the nodes next to it.
    # Contours is a dict, keyed by the node.id.
    cdef int i, j, off_i, off_j, oi, oj, ax, ax0, ax1, n1, n2
    cdef np.int64_t c1, c2
    cdef Node adj_node
    cdef VolumeContainer *vc1
    cdef VolumeContainer *vc0 = vcs[nid]
    cdef int s = (vc0.dims[1]*vc0.dims[0]
                + vc0.dims[0]*vc0.dims[2]
                + vc0.dims[1]*vc0.dims[2]) * 18
    # We allocate an array of fixed (maximum) size
    cdef np.ndarray[np.int64_t, ndim=2] joins = np.zeros((s, 2), dtype="int64")
    cdef int ti = 0, side, m1, m2, index
    cdef int pos[3]
    cdef int my_pos[3]
    cdef np.float64_t spos[3]

    for ax in range(3):
        ax0 = (ax + 1) % 3
        ax1 = (ax + 2) % 3
        n1 = vc0.dims[ax0]
        n2 = vc0.dims[ax1]
        for i in range(n1):
            for j in range(n2):
                for off_i in range(3):
                    oi = off_i - 1
                    if i == 0 and oi == -1: continue
                    if i == n1 - 1 and oi == 1: continue
                    for off_j in range(3):
                        oj = off_j - 1
                        if j == 0 and oj == -1: continue
                        if j == n2 - 1 and oj == 1: continue
                        pos[ax0] = i + oi
                        pos[ax1] = j + oj
                        my_pos[ax0] = i
                        my_pos[ax1] = j
                        for side in range(2):
                            # We go off each end of the block.
                            if side == 0:
                                pos[ax] = -1
                                my_pos[ax] = 0
                            else:
                                pos[ax] = vc0.dims[ax]
                                my_pos[ax] = vc0.dims[ax]-1
                            get_spos(vc0, pos[0], pos[1], pos[2], ax, spos)
                            adj_node = trunk._find_node(spos)
                            vc1 = vcs[adj_node.node_ind]
                            if spos_contained(vc1, spos):
                                index = vc_index(vc0, my_pos[0],
                                                 my_pos[1], my_pos[2])
                                m1 = vc0.mask[index]
                                c1 = (<np.int64_t*>vc0.data[0])[index]
                                index = vc_pos_index(vc1, spos)
                                m2 = vc1.mask[index]
                                c2 = (<np.int64_t*>vc1.data[0])[index]
                                if m1 == 1 and m2 == 1 and c1 > -1 and c2 > -1:
                                    if examined[adj_node.node_ind] == 0:
                                        joins[ti,0] = i64max(c1,c2)
                                        joins[ti,1] = i64min(c1,c2)
                                    else:
                                        joins[ti,0] = c1
                                        joins[ti,1] = c2
                                    ti += 1

    if ti == 0: return
    new_joins = tree.cull_joins(joins[:ti,:])
    tree.add_joins(new_joins)

@cython.boundscheck(False)
@cython.wraparound(False)
def update_joins(np.ndarray[np.int64_t, ndim=2] joins,
                 np.ndarray[np.int64_t, ndim=3] contour_ids,
                 np.ndarray[np.int64_t, ndim=1] final_joins):
    cdef int j, nj, nf
    cdef int ci, cj, ck
    nj = joins.shape[0]
    nf = final_joins.shape[0]
    for ci in range(contour_ids.shape[0]):
        for cj in range(contour_ids.shape[1]):
            for ck in range(contour_ids.shape[2]):
                if contour_ids[ci,cj,ck] == -1: continue
                for j in range(nj):
                    if contour_ids[ci,cj,ck] == joins[j,0]:
                        contour_ids[ci,cj,ck] = joins[j,1]
                        break
                for j in range(nf):
                    if contour_ids[ci,cj,ck] == final_joins[j]:
                        contour_ids[ci,cj,ck] = j + 1
                        break

cdef class FOFNode:
    cdef np.int64_t tag, count
    def __init__(self, np.int64_t tag):
        self.tag = tag
        self.count = 0

cdef class ParticleContourTree(ContourTree):
    cdef np.float64_t linking_length, linking_length2
    cdef np.float64_t DW[3]
    cdef np.float64_t DLE[3]
    cdef np.float64_t DRE[3]
    cdef bint periodicity[3]
    cdef int minimum_count

    def __init__(self, linking_length, periodicity = (True, True, True),
                 int minimum_count = 8):
        cdef int i
        self.linking_length = linking_length
        self.linking_length2 = linking_length * linking_length
        self.first = self.last = NULL
        for i in range(3):
            self.periodicity[i] = periodicity[i]
        self.minimum_count = minimum_count

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def identify_contours(self, OctreeContainer octree,
                                np.ndarray[np.int64_t, ndim=1] dom_ind,
                                np.ndarray[anyfloat, ndim=2] positions,
                                np.ndarray[np.int64_t, ndim=1] particle_ids,
                                int domain_id, int domain_offset):
        cdef np.ndarray[np.int64_t, ndim=1] pdoms, pcount, pind, doff
        cdef np.float64_t pos[3]
        cdef Oct *oct = NULL
        cdef Oct **neighbors = NULL
        cdef OctInfo oi
        cdef ContourID *c0
        cdef np.int64_t moff = octree.get_domain_offset(domain_id + domain_offset)
        cdef np.int64_t i, j, k, n, nneighbors = -1, pind0, offset
        cdef int counter = 0
        cdef int verbose = 0
        pcount = np.zeros_like(dom_ind)
        doff = np.zeros_like(dom_ind) - 1
        # First, we find the oct for each particle.
        pdoms = np.zeros(positions.shape[0], dtype="int64")
        pdoms -= -1
        # First we allocate our container
        cdef ContourID **container = <ContourID**> malloc(
            sizeof(ContourID*) * positions.shape[0])
        for i in range(3):
            self.DW[i] = (octree.DRE[i] - octree.DLE[i])
            self.DLE[i] = octree.DLE[i]
            self.DRE[i] = octree.DRE[i]
        for i in range(positions.shape[0]):
            counter += 1
            container[i] = NULL
            for j in range(3):
                pos[j] = positions[i, j]
            oct = octree.get(pos, NULL)
            if oct == NULL or (domain_id > 0 and oct.domain != domain_id):
                continue
            offset = oct.domain_ind - moff
            pcount[offset] += 1
            pdoms[i] = offset
        pind = np.argsort(pdoms)
        cdef np.int64_t *ipind = <np.int64_t*> pind.data
        cdef anyfloat *fpos = <anyfloat*> positions.data
        # pind is now the pointer into the position and particle_ids array.
        for i in range(positions.shape[0]):
            offset = pdoms[pind[i]]
            if doff[offset] < 0:
                doff[offset] = i
        del pdoms
        cdef int nsize = 27
        cdef np.int64_t *nind = <np.int64_t *> malloc(sizeof(np.int64_t)*nsize)
        counter = 0
        cdef np.int64_t frac = <np.int64_t> (doff.shape[0] / 20.0)
        if verbose == 1: print >> sys.stderr, "Will be outputting every", frac
        for i in range(doff.shape[0]):
            if verbose == 1 and counter >= frac:
                counter = 0
                print >> sys.stderr, "FOF-ing % 5.1f%% done" % ((100.0 * i)/doff.size)
            counter += 1
            # Any particles found for this oct?
            if doff[i] < 0: continue
            offset = pind[doff[i]]
            # This can probably be replaced at some point with a faster lookup.
            for j in range(3):
                pos[j] = positions[offset, j]
            oct = octree.get(pos, &oi)
            if oct == NULL or (domain_id > 0 and oct.domain != domain_id):
                continue
            # Now we have our primary oct, so we will get its neighbors.
            neighbors = octree.neighbors(&oi, &nneighbors, oct,
                                self.periodicity)
            # Now we have all our neighbors.  And, we should be set for what
            # else we need to do.
            if nneighbors > nsize:
                nind = <np.int64_t *> realloc(
                    nind, sizeof(np.int64_t)*nneighbors)
                nsize = nneighbors
            for j in range(nneighbors):
                nind[j] = neighbors[j].domain_ind - moff
                for n in range(j):
                    if nind[j] == nind[n]:
                        nind[j] = -1
                    break
            # This is allocated by the neighbors function, so we deallocate it.
            free(neighbors)
            # We might know that all our internal particles are linked.
            # Otherwise, we look at each particle.
            for j in range(pcount[i]):
                # Note that this offset is the particle index
                pind0 = pind[doff[i] + j]
                # Look at each neighboring oct
                for k in range(nneighbors):
                    if nind[k] == -1: continue
                    offset = doff[nind[k]]
                    if offset < 0: continue
                    # NOTE: doff[i] will not monotonically increase.  So we
                    # need a unique ID for each container that we are
                    # accessing.
                    self.link_particles(container,
                                        fpos, ipind,
                                        pcount[nind[k]],
                                        offset, pind0,
                                        doff[i] + j)
        cdef np.ndarray[np.int64_t, ndim=1] contour_ids
        contour_ids = np.ones(positions.shape[0], dtype="int64")
        contour_ids *= -1
        # Perform one last contour_find on each.  Note that we no longer need
        # to look at any of the doff or internal offset stuff.
        for i in range(positions.shape[0]):
            if container[i] == NULL: continue
            container[i] = contour_find(container[i])
        for i in range(positions.shape[0]):
            if container[i] == NULL: continue
            c0 = container[i]
            if c0.count < self.minimum_count: continue
            contour_ids[i] = particle_ids[pind[c0.contour_id]]
        free(container)
        del pind
        return contour_ids

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void link_particles(self, ContourID **container,
                                   anyfloat *positions,
                                   np.int64_t *pind,
                                   np.int64_t pcount,
                                   np.int64_t noffset,
                                   np.int64_t pind0,
                                   np.int64_t poffset):
        # Now we look at each particle and evaluate it
        cdef np.float64_t pos0[3]
        cdef np.float64_t pos1[3]
        cdef np.float64_t edges[2][3]
        cdef int link
        cdef ContourID *c0
        cdef ContourID *c1
        cdef np.int64_t pind1
        cdef int i, j
        # We use pid here so that we strictly take new ones.
        # Note that pind0 will not monotonically increase, but
        c0 = container[pind0]
        if c0 == NULL:
            c0 = container[pind0] = contour_create(poffset, self.last)
            self.last = c0
            if self.first == NULL:
                self.first = c0
        c0 = container[pind0] = contour_find(c0)
        for i in range(3):
            # We make a very conservative guess here about the edges.
            pos0[i] = positions[pind0*3 + i]
            edges[0][i] = pos0[i] - self.linking_length*1.01
            edges[1][i] = pos0[i] + self.linking_length*1.01
            if edges[0][i] < self.DLE[i] or edges[0][i] > self.DRE[i]:
                # We skip this one, since we're close to the boundary
                edges[0][i] = -1e30
                edges[1][i] = 1e30
        # Lets set up some bounds for the particles.  Maybe we can get away
        # with reducing our number of calls to r2dist_early.
        for i in range(pcount):
            pind1 = pind[noffset + i]
            if pind1 == pind0: continue
            c1 = container[pind1]
            if c1 != NULL and c1.contour_id == c0.contour_id:
                # Already linked.
                continue
            for j in range(3):
                pos1[j] = positions[pind1*3 + j]
            link = r2dist_early(pos0, pos1, self.DW, self.periodicity,
                                self.linking_length2, edges)
            if link == 0: continue
            if c1 == NULL:
                c0.count += 1
                container[pind1] = c0
            elif c0.contour_id != c1.contour_id:
                contour_union(c0, c1)
                c0 = container[pind1] = container[pind0] = contour_find(c0)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int r2dist_early(np.float64_t ppos[3],
                             np.float64_t cpos[3],
                             np.float64_t DW[3],
                             bint periodicity[3],
                             np.float64_t max_r2,
                             np.float64_t edges[2][3]):
    cdef int i
    cdef np.float64_t r2, DR
    r2 = 0.0
    for i in range(3):
        if cpos[i] < edges[0][i]:
            return 0
        if cpos[i] > edges[1][i]:
            return 0
    for i in range(3):
        DR = (ppos[i] - cpos[i])
        if not periodicity[i]:
            pass
        elif (DR > DW[i]/2.0):
            DR -= DW[i]
        elif (DR < -DW[i]/2.0):
            DR += DW[i]
        r2 += DR * DR
        if r2 > max_r2: return 0
    return 1
