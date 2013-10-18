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
from libc.stdlib cimport malloc, free

from amr_kdtools cimport _find_node, Node
from grid_traversal cimport VolumeContainer, PartitionedGrid, \
    vc_index, vc_pos_index

cdef extern from "math.h":
    double fabs(double x)

cdef extern from "stdlib.h":
    # NOTE that size_t might not be int
    void *alloca(int)

cdef inline np.int64_t i64max(np.int64_t i0, np.int64_t i1):
    if i0 > i1: return i0
    return i1

cdef inline np.int64_t i64min(np.int64_t i0, np.int64_t i1):
    if i0 < i1: return i0
    return i1

cdef struct ContourID

cdef struct ContourID:
    np.int64_t contour_id
    ContourID *parent
    ContourID *next
    ContourID *prev

cdef ContourID *contour_create(np.int64_t contour_id,
                               ContourID *prev = NULL):
    node = <ContourID *> malloc(sizeof(ContourID))
    #print "Creating contour with id", contour_id
    node.contour_id = contour_id
    node.next = node.parent = NULL
    node.prev = prev
    if prev != NULL: prev.next = node
    return node

cdef void contour_delete(ContourID *node):
    if node.prev != NULL: node.prev.next = node.next
    if node.next != NULL: node.next.prev = node.prev
    free(node)

cdef ContourID *contour_find(ContourID *node):
    cdef ContourID *temp, *root
    root = node
    while root.parent != NULL and root.parent != root:
        root = root.parent
    if root == root.parent: root.parent = NULL
    while node.parent != NULL:
        temp = node.parent
        node.parent = root
        node = temp
    return root

cdef void contour_union(ContourID *node1, ContourID *node2):
    if node1.contour_id < node2.contour_id:
        node2.parent = node1
    elif node2.contour_id < node1.contour_id:
        node1.parent = node2

cdef struct CandidateContour

cdef struct CandidateContour:
    np.int64_t contour_id
    np.int64_t join_id
    CandidateContour *next

cdef int candidate_contains(CandidateContour *first,
                            np.int64_t contour_id,
                            np.int64_t join_id = -1):
    while first != NULL:
        if first.contour_id == contour_id \
            and first.join_id == join_id: return 1
        first = first.next
    return 0

cdef CandidateContour *candidate_add(CandidateContour *first,
                                     np.int64_t contour_id,
                                     np.int64_t join_id = -1):
    cdef CandidateContour *node
    node = <CandidateContour *> malloc(sizeof(CandidateContour))
    node.contour_id = contour_id
    node.join_id = join_id
    node.next = first
    return node

cdef class ContourTree:
    cdef ContourID *first
    cdef ContourID *last

    def clear(self):
        # Here, we wipe out ALL of our contours, but not the pointers to them
        cdef ContourID *cur, *next
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
        # This is a helper function.
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
        while first != NULL:
            contours[i] = first.contour_id
            i += 1
            temp = first.next
            free(first)
            first = temp
        return contours

    def cull_joins(self, np.ndarray[np.int64_t, ndim=2] cjoins):
        cdef int i, j, k, ni, nj, nk, nc
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
        cdef ContourID *cur, *root
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
        cdef ContourID *cur, *root
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
                                np.int64_t start):
        cdef int i, j, k, ni, nj, nk, offset
        cdef int off_i, off_j, off_k, oi, ok, oj
        cdef ContourID *cur = NULL
        cdef ContourID *c1, *c2
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
                                if off_k > k and off_j > j and off_i > i:
                                    continue
                                off_k = ok - 1 + k
                                if not (0 <= off_k < nk): continue
                                offset = off_i*nj*nk + off_j*nk + off_k
                                c2 = container[offset]
                                if c2 == NULL: continue
                                c2 = contour_find(c2)
                                contour_union(cur, c2)
                                cur = contour_find(cur)
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    c1 = container[i*nj*nk + j*nk + k]
                    if c1 == NULL: continue
                    cur = c1
                    c1 = contour_find(c1)
                    contour_ids[i,j,k] = c1.contour_id
        
        for i in range(ni*nj*nk): 
            if container[i] != NULL: free(container[i])
        free(container)

#@cython.boundscheck(False)
@cython.wraparound(False)
def link_node_contours(Node trunk, contours, ContourTree tree):
    cdef int n_nodes = max(contours)
    cdef VolumeContainer **vcs = <VolumeContainer **> malloc(
        sizeof(VolumeContainer*) * n_nodes)
    cdef int i
    cdef PartitionedGrid pg
    for i in range(n_nodes):
        v = contours.get(i, None)
        if v is None:
            vcs[i] = NULL
            continue
        pg = v
        vcs[i] = pg.container
    cdef np.ndarray[np.uint8_t] examined = np.zeros(n_nodes, "uint8")
    for nid, (level, pg) in sorted(contours.items(), key = lambda a: -a[1][0]):
        construct_boundary_relationships(trunk, tree, nid, examined,
            vcs)

cdef inline void get_spos(VolumeContainer *vc, int i, int j, int k,
                          int axis, np.float64_t *spos):
    spos[0] = vc.left_edge[0] + i * vc.dds[0]
    spos[1] = vc.left_edge[1] + j * vc.dds[1]
    spos[2] = vc.left_edge[2] + k * vc.dds[2]
    spos[axis] += 0.5 * vc.dds[axis]

cdef inline int spos_contained(VolumeContainer *vc, np.float64_t *spos):
    cdef int i
    for i in range(3):
        if spos[i] < vc.left_edge[i] or spos[i] > vc.right_edge[i]: return 0
    return 1

#@cython.boundscheck(False)
@cython.wraparound(False)
cdef void construct_boundary_relationships(Node trunk, ContourTree tree, 
                np.int64_t nid, np.ndarray[np.uint8_t, ndim=1] examined,
                VolumeContainer **vcs):
    # We only look at the boundary and find the nodes next to it.
    # Contours is a dict, keyed by the node.id.
    cdef int i, j, nx, ny, nz, offset_i, offset_j, oi, oj, level
    cdef np.int64_t c1, c2
    cdef Node adj_node
    cdef VolumeContainer *vc1, *vc0 = vcs[nid]
    nx = vc0.dims[0]
    ny = vc0.dims[1]
    nz = vc0.dims[2]
    cdef int s = (ny*nx + nx*nz + ny*nz) * 18
    # We allocate an array of fixed (maximum) size
    cdef np.ndarray[np.int64_t, ndim=2] joins = np.zeros((s, 2), dtype="int64")
    cdef int ti = 0
    cdef int index
    cdef np.float64_t spos[3]

    # First the x-pass
    for i in range(ny):
        for j in range(nz):
            for offset_i in range(3):
                oi = offset_i - 1
                for offset_j in range(3):
                    oj = offset_j - 1
                    # Adjust by -1 in x, then oi and oj in y and z
                    get_spos(vc0, -1, i + oi, j + oj, 0, spos)
                    adj_node = _find_node(trunk, spos)
                    vc1 = vcs[adj_node.node_id]
                    if examined[adj_node.node_id] == 0 and \
                       spos_contained(vc1, spos):
                        # This is outside our VC, as 0 is a boundary layer
                        index = vc_index(vc0, 0, i, j)
                        c1 = (<np.int64_t*>vc0.data[0])[index]
                        index = vc_pos_index(vc1, spos)
                        c2 = (<np.int64_t*>vc1.data[0])[index]
                        if c1 > -1 and c2 > -1:
                            joins[ti,0] = i64max(c1,c2)
                            joins[ti,1] = i64min(c1,c2)
                            ti += 1
                    # This is outside our vc
                    get_spos(vc0, nx, i + oi, j + oj, 0, spos)
                    adj_node = _find_node(trunk, spos)
                    vc1 = vcs[adj_node.node_id]
                    if examined[adj_node.node_id] == 0 and \
                       spos_contained(vc1, spos):
                        # This is outside our VC, as 0 is a boundary layer
                        index = vc_index(vc0, nx - 1, i, j)
                        c1 = (<np.int64_t*>vc0.data[0])[index]
                        index = vc_pos_index(vc1, spos)
                        c2 = (<np.int64_t*>vc1.data[0])[index]
                        if c1 > -1 and c2 > -1:
                            joins[ti,0] = i64max(c1,c2)
                            joins[ti,1] = i64min(c1,c2)
                            ti += 1
    # Now y-pass
    for i in range(nx):
        for j in range(nz):
            for offset_i in range(3):
                oi = offset_i - 1
                if i == 0 and oi == -1: continue
                if i == nx - 1 and oi == 1: continue
                for offset_j in range(3):
                    oj = offset_j - 1
                    get_spos(vc0, i + oi, -1, j + oj, 1, spos)
                    adj_node = _find_node(trunk, spos)
                    vc1 = vcs[adj_node.node_id]
                    if examined[adj_node.node_id] == 0 and \
                       spos_contained(vc1, spos):
                        # This is outside our VC, as 0 is a boundary layer
                        index = vc_index(vc0, i, 0, j)
                        c1 = (<np.int64_t*>vc0.data[0])[index]
                        index = vc_pos_index(vc1, spos)
                        c2 = (<np.int64_t*>vc1.data[0])[index]
                        if c1 > -1 and c2 > -1:
                            joins[ti,0] = i64max(c1,c2)
                            joins[ti,1] = i64min(c1,c2)
                            ti += 1
                    get_spos(vc0, i + oi, ny, j + oj, 1, spos)
                    adj_node = _find_node(trunk, spos)
                    vc1 = vcs[adj_node.node_id]
                    if examined[adj_node.node_id] == 0 and \
                       spos_contained(vc1, spos):
                        # This is outside our VC, as 0 is a boundary layer
                        index = vc_index(vc0, i, ny, j)
                        c1 = (<np.int64_t*>vc0.data[0])[index]
                        index = vc_pos_index(vc1, spos)
                        c2 = (<np.int64_t*>vc1.data[0])[index]
                        if c1 > -1 and c2 > -1:
                            joins[ti,0] = i64max(c1,c2)
                            joins[ti,1] = i64min(c1,c2)
                            ti += 1

    # Now z-pass
    for i in range(nx):
        for j in range(ny):
            for offset_i in range(3):
                oi = offset_i - 1
                for offset_j in range(3):
                    oj = offset_j - 1
                    get_spos(vc0, i + oi,  j + oj, -1, 2, spos)
                    adj_node = _find_node(trunk, spos)
                    vc1 = vcs[adj_node.node_id]
                    if examined[adj_node.node_id] == 0 and \
                       spos_contained(vc1, spos):
                        # This is outside our VC, as 0 is a boundary layer
                        index = vc_index(vc0, i, j, 0)
                        c1 = (<np.int64_t*>vc0.data[0])[index]
                        index = vc_pos_index(vc1, spos)
                        c2 = (<np.int64_t*>vc1.data[0])[index]

                    get_spos(vc0, i + oi, j + oj, nz, 2, spos)
                    adj_node = _find_node(trunk, spos)
                    vc1 = vcs[adj_node.node_id]
                    if examined[adj_node.node_id] == 0 and \
                       spos_contained(vc1, spos):
                        # This is outside our VC, as 0 is a boundary layer
                        index = vc_index(vc0, i, j, nz)
                        c1 = (<np.int64_t*>vc0.data[0])[index]
                        index = vc_pos_index(vc1, spos)
                        c2 = (<np.int64_t*>vc1.data[0])[index]
                        if c1 > -1 and c2 > -1:
                            joins[ti,0] = i64max(c1,c2)
                            joins[ti,1] = i64min(c1,c2)
                            ti += 1
    new_joins = tree.cull_joins(tree[:ti,:])
    tree.add_joins(new_joins)

cdef inline int are_neighbors(
            np.float64_t x1, np.float64_t y1, np.float64_t z1,
            np.float64_t dx1, np.float64_t dy1, np.float64_t dz1,
            np.float64_t x2, np.float64_t y2, np.float64_t z2,
            np.float64_t dx2, np.float64_t dy2, np.float64_t dz2,
        ):
    # We assume an epsilon of 1e-15
    if fabs(x1-x2) > 0.5*(dx1+dx2): return 0
    if fabs(y1-y2) > 0.5*(dy1+dy2): return 0
    if fabs(z1-z2) > 0.5*(dz1+dz2): return 0
    return 1

@cython.boundscheck(False)
@cython.wraparound(False)
def identify_field_neighbors(
            np.ndarray[dtype=np.float64_t, ndim=1] field,
            np.ndarray[dtype=np.float64_t, ndim=1] x,
            np.ndarray[dtype=np.float64_t, ndim=1] y,
            np.ndarray[dtype=np.float64_t, ndim=1] z,
            np.ndarray[dtype=np.float64_t, ndim=1] dx,
            np.ndarray[dtype=np.float64_t, ndim=1] dy,
            np.ndarray[dtype=np.float64_t, ndim=1] dz,
        ):
    # We assume this field is pre-jittered; it has no identical values.
    cdef int outer, inner, N, added
    cdef np.float64_t x1, y1, z1, dx1, dy1, dz1
    N = field.shape[0]
    #cdef np.ndarray[dtype=np.object_t] joins
    joins = [[] for outer in range(N)]
    #joins = np.empty(N, dtype='object')
    for outer in range(N):
        if (outer % 10000) == 0: print outer, N
        x1 = x[outer]
        y1 = y[outer]
        z1 = z[outer]
        dx1 = dx[outer]
        dy1 = dy[outer]
        dz1 = dz[outer]
        this_joins = joins[outer]
        added = 0
        # Go in reverse order
        for inner in range(outer, 0, -1):
            if not are_neighbors(x1, y1, z1, dx1, dy1, dz1,
                                 x[inner], y[inner], z[inner],
                                 dx[inner], dy[inner], dz[inner]):
                continue
            # Hot dog, we have a weiner!
            this_joins.append(inner)
            added += 1
            if added == 26: break
    return joins

@cython.boundscheck(False)
@cython.wraparound(False)
def extract_identified_contours(int max_ind, joins):
    cdef int i
    contours = []
    for i in range(max_ind + 1): # +1 to get to the max_ind itself
        contours.append(set([i]))
        if len(joins[i]) == 0:
            continue
        proto_contour = [i]
        for j in joins[i]:
            proto_contour += contours[j]
        proto_contour = set(proto_contour)
        for j in proto_contour:
            contours[j] = proto_contour
    return contours

@cython.boundscheck(False)
@cython.wraparound(False)
def update_joins(np.ndarray[np.int64_t, ndim=2] joins,
                 np.ndarray[np.int64_t, ndim=1] contour_ids):
    cdef np.int64_t new, old
    cdef int i, j, nc, nj
    nc = contour_ids.shape[0]
    nj = joins.shape[0]
    for i in range(nc):
        if contour_ids[i] == -1: continue
        for j in range(nj):
            if contour_ids[i] == joins[j,0]:
                contour_ids[i] = joins[j,1]
                break
