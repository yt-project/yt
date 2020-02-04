cimport numpy as numpy
import numpy as np
from libcpp cimport bool
from libc.math cimport floor
cimport cython

from yt.geometry.oct_container cimport SparseOctreeContainer, OctInfo
from yt.geometry.oct_visitors cimport Oct, IndexOcts, StoreIndex
from yt.geometry.selection_routines cimport SelectorObject, AlwaysSelector


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
    cdef np.float64_t tx0, tx1, ty0, ty1, tz0, tz1, tmin, tmax
    cdef Oct *oct
    cdef int i
    cdef int ind[3]
    cdef np.float64_t dds[3], 
    cdef list octList, cellList
    oct = NULL

    for i in range(3):
        o[i] = 0  # origin
        d[i] = 0  # direction
    
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

    print('a=%s' % a)
    # Compute interesections with all 6 planes of node
    tx0 = (octree.DLE[0] - o[0]) / d[0]
    tx1 = (octree.DRE[0] - o[0]) / d[0]
    ty0 = (octree.DLE[1] - o[1]) / d[1]
    ty1 = (octree.DRE[1] - o[1]) / d[1]
    tz0 = (octree.DLE[2] - o[2]) / d[2]
    tz1 = (octree.DRE[2] - o[2]) / d[2]

    tmin = max(tx0, ty0, tz0)
    tmax = min(tx1, ty1, tz1)

    rr = r.at(tmin)
    # TODO: call for all roots
    ii = 1
    for i in range(3):
        dds[i] = (octree.DRE[i] - octree.DLE[i])/octree.nn[i]
        ind[i] = <int> (floor((rr[i] - octree.DLE[i])/dds[i]))
        if a & ii:
            ind[i] = octree.nn[i] - ind[i]
        ii <<= 1
    octree.get_root(ind, &oct)

    if oct == NULL:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)

    # Hits, so process subtree
    octList = []
    cellList = []
    if (tmin < tmax) and (tmax > 0):
        proc_subtree(tx0, ty0, tz0, tx1, ty1, tz1, oct, a, octList, cellList)
    return np.asarray(octList, dtype=np.int64), np.asarray(cellList, dtype=np.int64)

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
    # Find entry plane
    tmax = max(tx0, ty0, tz0)
    if tmax == tx0:
        entry_plane = YZ
    elif tmax == ty0:
        entry_plane = XZ
    elif tmax == tz0:
        entry_plane = XY
    cdef np.uint8_t first_node

    # Now find first node
    first_node = 0
    if entry_plane == XY:
        if txM < tz0:
            first_node |= 0b001
        if tyM < tz0:
            first_node |= 0b010
    elif entry_plane == YZ:
        if tyM < tx0:
            first_node |= 0b010
        if tzM < tx0:
            first_node |= 0b100
    elif entry_plane == XZ:
        if txM < ty0:
            first_node |= 0b001
        if tzM < ty0:
            first_node |= 0b100
    return first_node

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.uint8_t next_Node(np.uint8_t currNode,
        np.float64_t txM, np.float64_t tyM, np.float64_t tzM,
        np.float64_t tx1, np.float64_t ty1, np.float64_t tz1):
    cdef np.float64_t tmin, tx, ty, tz

    if currNode & 0b001:
        tx = tx1
    else:
        tx = txM

    if currNode & 0b010:
        ty = ty1
    else:
        ty = tyM

    if currNode & 0b100:
        tz = tz1
    else:
        tz = tzM

    tmin = min(tx, ty, tz)

    if tmin == tx:    # Next node in x direction
        if currNode & 0b001:
            return 8
        else:
            return currNode + 1
    elif tmin == ty:  # Next node in y direction
        if currNode & 0b010:
            return 8
        else:
            return currNode + 2
    else:             # Next node in z direction
        if currNode & 0b100:
            return 8
        else:
            return currNode + 4

cdef inline bool isLeaf(const Oct *o):
    return o.children == NULL

# TODO: support negative directions
# TODO: support ray length
cdef void proc_subtree(
        np.float64_t tx0, np.float64_t ty0, np.float64_t tz0,
        np.float64_t tx1, np.float64_t ty1, np.float64_t tz1,
        Oct* oct, int a, list octList, list cellList, int level=0):
    cdef np.uint8_t currNode
    cdef np.float64_t txM, tyM, tzM
    cdef np.uint8_t entry_plane, exit_plane
    cdef bool leaf
    if tx1 < 0 or ty1 < 0 or tz1 < 0:
        return

    leaf = isLeaf(oct)

    # Compute midpoints
    txM = (tx0 + tx1) / 2.
    tyM = (ty0 + ty1) / 2.
    tzM = (tz0 + tz1) / 2.

    # Compute entry/exit planes
    # entry_plane = find_entry_plane(tx0, ty0, tz0)

    currNode = find_firstNode(tx0, ty0, tz0, txM, tyM, tzM)

    # print('%s[%s]@lvl=%s' % ('\t'*level, oct.domain_ind, currNode, level))

    while True:
        # print('%scurrNode=%s' %('\t'*level, currNode))
        # print('%scurrNode=%s %s %.2f %.2f %.2f %.2f %.2f %.2f' % ('\t'*level, currNode, a, txM, tyM, tzM, tx1, ty1, tz1))

        if leaf:
            octList.append(oct.domain_ind)
            cellList.append(currNode^a)
        else:
            if currNode == 0:
                proc_subtree(tx0, ty0, tz0, txM, tyM, tzM, oct.children[a], a, octList, cellList, level+1)
            elif currNode == 1:
                proc_subtree(txM, ty0, tz0, tx1, tyM, tzM, oct.children[4^a], a, octList, cellList, level+1)
            elif currNode == 2:
                proc_subtree(tx0, tyM, tz0, txM, ty1, tzM, oct.children[2^a], a, octList, cellList, level+1)
            elif currNode == 3:
                proc_subtree(txM, tyM, tz0, tx1, ty1, tzM, oct.children[6^a], a, octList, cellList, level+1)
            elif currNode == 4:
                proc_subtree(tx0, ty0, tzM, txM, tyM, tz1, oct.children[1^a], a, octList, cellList, level+1)
            elif currNode == 5:
                proc_subtree(txM, ty0, tzM, tx1, tyM, tz1, oct.children[5^a], a, octList, cellList, level+1)
            elif currNode == 6:
                proc_subtree(tx0, tyM, tzM, txM, ty1, tz1, oct.children[3^a], a, octList, cellList, level+1)
            elif currNode == 7:
                proc_subtree(txM, tyM, tzM, tx1, ty1, tz1, oct.children[7], a, octList, cellList, level+1)
            
        currNode = next_Node(currNode, txM, tyM, tzM, tx1, ty1, tz1)

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

    return cell_inds