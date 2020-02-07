"""
Adapted from "An Efﬁcient Parametric Algorithm for Octree Traversal", J. Revelles, C. Ureña, M. Lastra, 2000
http://wscg.zcu.cz/wscg2000/Papers_2000/X31.pdf
"""
# distutils: language = c++

cimport numpy as np
import numpy as np
from libcpp cimport bool
from libcpp.vector cimport vector
from libc.math cimport floor
cimport cython
from cpython cimport Py_buffer

from yt.geometry.oct_container cimport SparseOctreeContainer, OctInfo
from yt.geometry.oct_visitors cimport Oct, IndexOcts, StoreIndex
from yt.geometry.selection_routines cimport SelectorObject, AlwaysSelector


cdef class Uint8VectorHolder:
    # See https://cython.readthedocs.io/en/latest/src/userguide/buffer.html#a-matrix-class
    cdef vector[np.uint8_t] v
    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]
    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Py_ssize_t itemsize = sizeof(self.v[0])

        self.shape[0] = self.v.size()

        # Stride 1 is the distance, in bytes, between two items in a row;
        # this is the distance between two adjacent items in the vector.
        # Stride 0 is the distance between the first elements of adjacent rows.
        self.strides[0] = <Py_ssize_t>(  <char *>&(self.v[1])
                                       - <char *>&(self.v[0]))

        buffer.buf = <char *>&(self.v[0])
        buffer.format = 'B'                     # unsigned long
        buffer.internal = NULL                  # see References
        buffer.itemsize = itemsize
        buffer.len = self.v.size() * itemsize   # product(shape) * itemsize
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

cdef class Uint64VectorHolder:
    # See https://cython.readthedocs.io/en/latest/src/userguide/buffer.html#a-matrix-class
    cdef vector[np.uint64_t] v
    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]
    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Py_ssize_t itemsize = sizeof(self.v[0])

        self.shape[0] = self.v.size()

        # Stride 1 is the distance, in bytes, between two items in a row;
        # this is the distance between two adjacent items in the vector.
        # Stride 0 is the distance between the first elements of adjacent rows.
        self.strides[0] = <Py_ssize_t>(  <char *>&(self.v[1])
                                       - <char *>&(self.v[0]))

        buffer.buf = <char *>&(self.v[0])
        buffer.format = 'L'                     # unsigned long
        buffer.internal = NULL                  # see References
        buffer.itemsize = itemsize
        buffer.len = self.v.size() * itemsize   # product(shape) * itemsize
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

cdef class Float64VectorHolder:
    # See https://cython.readthedocs.io/en/latest/src/userguide/buffer.html#a-matrix-class
    cdef vector[np.float64_t] v
    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]
    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Py_ssize_t itemsize = sizeof(self.v[0])

        self.shape[0] = self.v.size()

        # Stride 1 is the distance, in bytes, between two items in a row;
        # this is the distance between two adjacent items in the vector.
        # Stride 0 is the distance between the first elements of adjacent rows.
        self.strides[0] = <Py_ssize_t>(  <char *>&(self.v[1])
                                       - <char *>&(self.v[0]))

        buffer.buf = <char *>&(self.v[0])
        buffer.format = 'd'                     # float
        buffer.internal = NULL                  # see References
        buffer.itemsize = itemsize
        buffer.len = self.v.size() * itemsize   # product(shape) * itemsize
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

cdef class Ray(object):
    def __init__(self, np.ndarray origin, np.ndarray direction, np.float64_t length):
        self.origin = np.asarray(origin)
        self.direction = direction / np.linalg.norm(direction)
        self.length = length

    @property
    def end(self):
        return self.at(self.length)

    cpdef np.ndarray[np.float64_t, ndim=1] at(self, np.float64_t t):
        cdef np.ndarray[np.float64_t, ndim=1] out = np.empty(3)
        for i in range(3):
            out[i] = self.origin[i] + t * self.direction[i]
        return out

def ray_step(SparseOctreeContainer octree, Ray r):
    cdef np.uint8_t a = 0
    cdef np.int8_t ii
    cdef np.float64_t o[3]
    cdef np.float64_t rr[3]
    cdef np.float64_t d[3]
    cdef np.float64_t tx0, tx1, ty0, ty1, tz0, tz1, tmin_domain, tmax_domain
    cdef Oct *oct
    cdef int i
    cdef int ind[3]
    cdef np.float64_t dds[3]
    oct = NULL

    ii = 1
    for i in range(3):
        if r.direction[i] < 0:
            o[i] = octree.DRE[i] - r.origin[i]
            d[i] = max(1e-99, -r.direction[i])
            a |= ii
        else:
            o[i] = r.origin[i]
            d[i] = max(1e-99, r.direction[i])
        ii <<= 1

    # Compute intersections with all 6 planes of octree
    tx0 = (octree.DLE[0] - o[0]) / d[0]
    ty0 = (octree.DLE[1] - o[1]) / d[1]
    tz0 = (octree.DLE[2] - o[2]) / d[2]

    tx1 = (octree.DRE[0] - o[0]) / d[0]
    ty1 = (octree.DRE[1] - o[1]) / d[1]
    tz1 = (octree.DRE[2] - o[2]) / d[2]

    tmin_domain = max(tx0, ty0, tz0)
    tmax_domain = min(tx1, ty1, tz1)

    # No hit at all, early break
    if (tmin_domain > tmax_domain) or (tmax_domain < 0):
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64), np.array([], dtype=np.float64)

    # Containers
    cdef Uint64VectorHolder octList = Uint64VectorHolder()
    cdef Uint8VectorHolder cellList = Uint8VectorHolder()
    cdef Float64VectorHolder tList = Float64VectorHolder()

    cdef np.float64_t txin, tyin, tzin, dtx, dty, dtz, tmin, tmax, txout, tyout, tzout

    # Locate first node
    # FIXME: add small epsilon ?
    rr = r.at(tmin_domain)

    for i in range(3):
        dds[i] = (octree.DRE[i] - octree.DLE[i])/octree.nn[i]
        ind[i] = <int> (floor((rr[i] - octree.DLE[i])/dds[i]))

    # Compute local in/out t
    dtx = dds[0]/d[0]
    dty = dds[1]/d[1]
    dtz = dds[2]/d[2]

    txin = tx0+ind[0]*dtx
    tyin = ty0+ind[1]*dty
    tzin = tz0+ind[2]*dtz

    txout = txin+dtx
    tyout = tyin+dty
    tzout = tzin+dtz

    tmin = max(txin, tyin, tzin)
    tmax = min(txout, tyout, tzout)

    # Loop over all cells until reaching the out face
    while tmax < tmax_domain:
        octree.get_root(ind, &oct)

        if oct != NULL and tmax > 0:  # no need to check tmin < tmax, as dtx,dty,dtz > 0
            # Hits, so process subtree
            proc_subtree(txin, tyin, tzin, txout, tyout, tzout,
                oct, a, octList.v, cellList.v, tList.v)
        if tmax == txout:    # Next x
            ind[0] += 1
            txin = txout
            txout += dtx
        elif tmax == tyout:  # Next y
            ind[1] += 1
            tyin = tyout
            tyout += dty
        elif tmax == tzout:  # Next z
            ind[2] += 1
            tzin = tzout
            tzout += dtz

        # Local in/out ts
        tmin = max(txin, tyin, tzin)
        tmax = min(txout, tyout, tzout)

    return np.asarray(octList), np.asarray(cellList), np.asarray(tList)

cdef np.uint8_t YZ = 0
cdef np.uint8_t XZ = 1
cdef np.uint8_t XY = 2

cdef np.uint8_t find_entry_plane(np.float64_t tx0, np.float64_t ty0, np.float64_t tz0):
    cdef np.float64_t tmax
    # Find entry plane
    tmax = max(tx0, ty0, tz0)
    if tmax == tx0:
        return YZ
    elif tmax == ty0:
        return XZ
    elif tmax == tz0:
        return XY


cdef np.uint8_t find_firstNode(
        np.float64_t tx0, np.float64_t ty0, np.float64_t tz0,
        np.float64_t txM,np.float64_t tyM, np.float64_t tzM):

    cdef np.float64_t tmax
    cdef np.uint8_t entry_plane
    cdef np.uint8_t first_node
    first_node = 0

    tmax = max(tx0, ty0, tz0)
    if tmax == tx0:    # YZ plane
        if tyM < tx0:
            first_node |= 0b010
        if tzM < tx0:
            first_node |= 0b001
    elif tmax == ty0:  # XZ plane
        if txM < ty0:
            first_node |= 0b100
        if tzM < ty0:
            first_node |= 0b001
    elif tmax == tz0:  # XY plane
        if txM < tz0:
            first_node |= 0b100
        if tyM < tz0:
            first_node |= 0b010

    return first_node

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.uint8_t next_node(np.uint8_t currNode,
        np.float64_t tx1, np.float64_t ty1, np.float64_t tz1):
    cdef np.float64_t tmin
    tmin = min(tx1, ty1, tz1)

    if tmin == tx1:    # YZ plane, increase x
        if currNode & 0b100:
            return 8
        else:
            return currNode + 4
    elif tmin == ty1:  # XZ plane, increase y
        if currNode & 0b010:
            return 8
        else:
            return currNode + 2
    else:             # XY plane, increase z
        if currNode & 0b001:
            return 8
        else:
            return currNode + 1

cdef inline bool isLeaf(const Oct *o, np.uint8_t currNode):
    return (o.children == NULL) or (o.children[currNode] == NULL)

cdef inline np.uint8_t swap3bits(const np.uint8_t lev):
    # Swap little-endian to big-endian (C to F ordering)
    return ((lev & 0b100) >> 2) + (lev & 0b010) + ((lev & 0b001) << 2)

# TODO: support negative directions
# TODO: support ray length
@cython.cdivision(True)
cdef void proc_subtree(
        const np.float64_t tx0, const np.float64_t ty0, const np.float64_t tz0,
        const np.float64_t tx1, const np.float64_t ty1, const np.float64_t tz1,
        const Oct* oct, const int a, vector[np.uint64_t] &octList, vector[np.uint8_t] &cellList, vector[np.float64_t] &tList,
        int level=0):
    cdef np.uint8_t currNode, nextNode
    cdef np.float64_t txM, tyM, tzM
    cdef np.uint8_t entry_plane, exit_plane
    cdef bool leaf
    if tx1 < 0 or ty1 < 0 or tz1 < 0:
        return

    if oct == NULL:
        print('This should not happen!')
        return

    # Compute midpoints
    txM = (tx0 + tx1) / 2.
    tyM = (ty0 + ty1) / 2.
    tzM = (tz0 + tz1) / 2.

    currNode = find_firstNode(tx0, ty0, tz0, txM, tyM, tzM)

    while True:
        leaf = isLeaf(oct, currNode^a)
        # print('%scurrNode=%s' %('\t'*level, currNode))
        # print('%scurrNode=%s %s (%.2f %.2f %.2f) (%.2f %.2f %.2f)' % ('\t'*level, currNode, a, txM, tyM, tzM, tx1, ty1, tz1))

        if leaf:
            octList.push_back(oct.domain_ind)
            # Need to swap bits before storing as octree is C-style in memory and F-style on file
            cellList.push_back(swap3bits(currNode^a))

        if currNode == 0:
            if not leaf: proc_subtree(tx0, ty0, tz0, txM, tyM, tzM, oct.children[  a], a, octList, cellList, tList, level+1)
            else: tList.push_back(min(txM, tyM, tzM))
            nextNode = next_node(currNode, txM, tyM, tzM)
        elif currNode == 1:
            if not leaf: proc_subtree(tx0, ty0, tzM, txM, tyM, tz1, oct.children[1^a], a, octList, cellList, tList, level+1)
            else: tList.push_back(min(txM, tyM, tz1))
            nextNode = next_node(currNode, txM, tyM, tz1)
        elif currNode == 2:
            if not leaf: proc_subtree(tx0, tyM, tz0, txM, ty1, tzM, oct.children[2^a], a, octList, cellList, tList, level+1)
            else: tList.push_back(min(txM, ty1, tzM))
            nextNode = next_node(currNode, txM, ty1, tzM)
        elif currNode == 3:
            if not leaf: proc_subtree(tx0, tyM, tzM, txM, ty1, tz1, oct.children[3^a], a, octList, cellList, tList, level+1)
            else: tList.push_back(min(txM, ty1, tz1))
            nextNode = next_node(currNode, txM, ty1, tz1)
        elif currNode == 4:
            if not leaf: proc_subtree(txM, ty0, tz0, tx1, tyM, tzM, oct.children[4^a], a, octList, cellList, tList, level+1)
            else: tList.push_back(min(tx1, tyM, tzM))
            nextNode = next_node(currNode, tx1, tyM, tzM)
        elif currNode == 5:
            if not leaf: proc_subtree(txM, ty0, tzM, tx1, tyM, tz1, oct.children[5^a], a, octList, cellList, tList, level+1)
            else: tList.push_back(min(tx1, tyM, tz1))
            nextNode = next_node(currNode, tx1, tyM, tz1)
        elif currNode == 6:
            if not leaf: proc_subtree(txM, tyM, tz0, tx1, ty1, tzM, oct.children[6^a], a, octList, cellList, tList, level+1)
            else: tList.push_back(min(tx1, ty1, tzM))
            nextNode = next_node(currNode, tx1, ty1, tzM)
        elif currNode == 7:
            if not leaf: proc_subtree(txM, tyM, tzM, tx1, ty1, tz1, oct.children[7  ], a, octList, cellList, tList, level+1)
            else: tList.push_back(min(tx1, ty1, tz1))
            nextNode = 8

        currNode = nextNode

        # Break when hitting 8'th node
        if currNode == 8:
            break

cpdef domain2ind(SparseOctreeContainer octree, SelectorObject selector,
            int domain_id=-1):
    cdef StoreIndex visitor
    cdef np.ndarray[np.int64_t, ndim=4] cell_inds

    cell_inds = np.empty((octree.nocts, 2, 2, 2), dtype=np.int64)

    visitor = StoreIndex(octree, domain_id)
    visitor.cell_inds = cell_inds

    octree.visit_all_octs(selector, visitor)

    return cell_inds.reshape(-1, 8)

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef octcells2ind(np.ndarray[np.int64_t, ndim=2] domain2ind, np.ndarray[np.uint64_t, ndim=1] oct_inds, np.ndarray[np.uint8_t, ndim=1] cell_inds):
    cdef np.ndarray[np.int64_t, ndim=1] indexes

    cdef np.int64_t[:] indexes_view
    cdef np.uint64_t[:] oct_inds_view = oct_inds
    cdef np.uint8_t[:] cell_inds_view = cell_inds
    cdef int i
    indexes = np.empty(oct_inds.size, dtype=np.int64)

    indexes_view = indexes

    for i in range(oct_inds.size):
        indexes_view[i] = domain2ind[oct_inds_view[i], cell_inds_view[i]]

    return indexes