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
from yt.frontends.ramses.io_utils cimport hilbert3d_single

cdef np.float64_t epsilon = 1.e-10  # TODO: find better size of epsilon. This corresponds to levelmax=33



cdef class Uint8VectorHolder:
    # See https://cython.readthedocs.io/en/latest/src/userguide/buffer.html#a-matrix-class
    cdef vector[np.uint8_t] v
    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]
    def __cinit__(self, int size):
        self.v.reserve(size)

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

cdef class Uint64VectorHolder:
    # See https://cython.readthedocs.io/en/latest/src/userguide/buffer.html#a-matrix-class
    cdef vector[np.uint64_t] v
    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]
    def __cinit__(self, int size):
        self.v.reserve(size)

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

cdef class Float64VectorHolder:
    # See https://cython.readthedocs.io/en/latest/src/userguide/buffer.html#a-matrix-class
    cdef vector[np.float64_t] v
    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]
    def __cinit__(self, int size):
        self.v.reserve(size)

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

cdef class Ray(object):
    def __init__(self, np.ndarray origin, np.ndarray direction, np.float64_t length):
        self.origin = np.asarray(origin)
        self.direction = direction / np.linalg.norm(direction)
        self.length = length

    @property
    def end(self):
        return self.at(self.length)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef np.ndarray[np.float64_t, ndim=1] at(self, np.float64_t t):
        cdef np.ndarray[np.float64_t, ndim=1] out = np.empty(3)
        for i in range(3):
            out[i] = self.origin[i] + t * self.direction[i]
        return out

    @cython.boundscheck(False) # turn of bounds-checking for entire function
    @cython.wraparound(False)  # turn of bounds-checking for entire function
    cpdef void trilinear(self, const np.float64_t tmin, const np.float64_t tmax,
                         const np.float64_t[:, :, :] data_in,
                         np.float64_t[:] data_out,
                         const np.float64_t[:] DLE,
                         const np.float64_t Deltax,
                         const int npt):
        """Interpolate npoint between tin and tout, given vertex-centred data.

        Parameters
        ----------
        tmin, tmax : float
        DLE : array (3,)
            The origin of the cell
        Deltax : float
            The size of the cell
        data_in : array (2, 2, 2)
            Vertex-centred data around the cell of interest
        data_out : array (npt, )
            The data (modified in-place)
        """
        cdef int i
        cdef np.float64_t dt
        cdef np.float64_t x0, y0, z0, x1, y1, z1, dx, dy, dz

        dt = (tmax-tmin)/npt

        dx = self.direction[0]/Deltax
        dy = self.direction[1]/Deltax
        dz = self.direction[2]/Deltax

        x0 = (self.origin[0]-DLE[0]) / Deltax + tmin*dx
        y0 = (self.origin[1]-DLE[1]) / Deltax + tmin*dy
        z0 = (self.origin[2]-DLE[2]) / Deltax + tmin*dz

        x1 = 1-x0
        y1 = 1-y0
        z1 = 1-z0

        for i in range(npt):
            # Tri-linear interpolation
            data_out[i] = (
                data_in[0,0,0] * x1 * y1 * z1 +
                data_in[1,0,0] * x0 * y1 * z1 +
                data_in[0,1,0] * x1 * y0 * z1 +
                data_in[1,1,0] * x0 * y0 * z1 +
                data_in[0,0,1] * x1 * y1 * z0 +
                data_in[1,0,1] * x0 * y1 * z0 +
                data_in[0,1,1] * x1 * y0 * z0 +
                data_in[1,1,1] * x0 * y0 * z0
            )

            x0 += dt*dx
            y0 += dt*dy
            z0 += dt*dz

            x1 = 1-x0
            y1 = 1-y0
            z1 = 1-z0

cdef class DomainFinder:
    cdef np.uint64_t[:] keys
    cdef np.float64_t bscale
    cdef int ncpu
    cdef int bit_length

    def __cinit__(self, ds):
        self.keys = ds.hilbert['keys'].astype(np.uint64)
        self.bscale = ds.hilbert['bscale']
        self.ncpu = ds.parameters['ncpu']
        self.bit_length = ds.hilbert['bit_length']

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int find_domain(self, np.float64_t[:] pos):
        cdef int nextDom
        cdef np.uint64_t ihilbert
        cdef np.int64_t ix, iy, iz

        ix = <np.uint64_t> (pos[0]*self.bscale)
        iy = <np.uint64_t> (pos[1]*self.bscale)
        iz = <np.uint64_t> (pos[2]*self.bscale)
        ihilbert = hilbert3d_single(ix, iy, iz, self.bit_length)

        # after the loop, nextDom contains the id of the first domain
        for nextDom in range(1, self.ncpu+1):
            if ihilbert < self.keys[nextDom]:
                break
        return nextDom


cpdef ray_step_multioctrees(dict octrees, Ray r, ds):
    # Find entry sparse octree
    cdef SparseOctreeContainer octree
    cdef int count, nextDom

    # Cell info containers -- preallocate size of domain
    octree = octrees[next(iter(octrees))]
    count = <int> (octree.nn[0] * np.sqrt(3))
    cdef Uint64VectorHolder octList = Uint64VectorHolder(count)
    cdef Uint8VectorHolder cellList = Uint8VectorHolder(count)
    cdef Float64VectorHolder tList = Float64VectorHolder(count)
    # Domain info containers -- preallocate number of domains
    count = <int> len(octrees)
    cdef Uint64VectorHolder countPerDomain = Uint64VectorHolder(count)
    cdef Uint64VectorHolder domainList = Uint64VectorHolder(count)

    # Find first domain using hilbert curve
    cdef int i
    cdef np.float64_t tmin, tin
    cdef np.float64_t[:] pos
    cdef DomainFinder df = DomainFinder(ds)

    tin = 0
    for i in range(3):
        if r.direction[i] < 0:
            tin = max(tin, octree.DRE[i] - r.origin[i]) / r.direction[i]
        else:
            tin = max(tin, octree.DLE[i] - r.origin[i]) / r.direction[i]

    pos = r.at(tin)

    nextDom = df.find_domain(pos)

    # Other variables
    cdef int nAdded
    tmin = 0
    count = 0
    nAdded = 0
    while nextDom > 0:
        # Add domain to list of domains
        domainList.v.push_back(nextDom)

        # Call ray traversal on domain
        octree = octrees[nextDom]
        nextDom = ray_step(octree, r, octList.v, cellList.v, tList.v, tmin, nextDom)

        # Update number of cells crossed
        nAdded = tList.v.size() - count
        count = tList.v.size()

        # Store this number
        countPerDomain.v.push_back(nAdded)

        # Next starting point
        if tList.v.size() > 0:
            tmin = tList.v.back()

    return np.asarray(octList), np.asarray(cellList), np.asarray(tList), np.asarray(domainList), np.asarray(countPerDomain)


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef int ray_step(SparseOctreeContainer octree, Ray r,
                   vector[np.uint64_t] &octList, vector[np.uint8_t] &cellList, vector[np.float64_t] &tList,
                   np.float64_t t0, const int curDom):
    cdef np.uint8_t a
    cdef np.int8_t ii
    cdef np.float64_t[3] o, rr, d
    cdef np.float64_t tx0, tx1, ty0, ty1, tz0, tz1, tmin_domain, tmax_domain
    cdef Oct *oct
    cdef int i, nextDom
    cdef int[3] ind, idx, aind
    cdef np.float64_t[3] dds

    nextDom = curDom
    oct = NULL
    ii = 0b100
    a = 0b000
    for i in range(3):
        if r.direction[i] < 0:
            o[i] = octree.DRE[i] - r.origin[i]
            d[i] = max(1e-99, -r.direction[i])
            a |= ii
            idx[i] = -1
            aind[i] = octree.nn[i]-1
        else:
            o[i] = r.origin[i]
            d[i] = max(1e-99, r.direction[i])
            idx[i] = +1
            aind[i] = 0
        ii >>= 1

    # Compute intersections with all 6 planes of octree
    tx0 = (octree.DLE[0] - o[0]) / d[0]
    ty0 = (octree.DLE[1] - o[1]) / d[1]
    tz0 = (octree.DLE[2] - o[2]) / d[2]

    tx1 = (octree.DRE[0] - o[0]) / d[0]
    ty1 = (octree.DRE[1] - o[1]) / d[1]
    tz1 = (octree.DRE[2] - o[2]) / d[2]

    tmin_domain = max(tx0, ty0, tz0, t0)
    tmax_domain = min(tx1, ty1, tz1)

    # No hit at all, early break
    if (tmin_domain > tmax_domain) or (tmax_domain < 0):
        return -1

    cdef np.float64_t txin, tyin, tzin, dtx, dty, dtz, tmin, tmax, txout, tyout, tzout

    # Locate first node
    rr = r.at(tmin_domain+epsilon)

    for i in range(3):
        dds[i] = (octree.DRE[i] - octree.DLE[i])/octree.nn[i]
        ind[i] = <int> (floor((rr[i] - octree.DLE[i])/dds[i]))

    # Compute local in/out t
    dtx = dds[0]/d[0]
    dty = dds[1]/d[1]
    dtz = dds[2]/d[2]

    txin = tx0+(ind[0]^aind[0])*dtx
    tyin = ty0+(ind[1]^aind[1])*dty
    tzin = tz0+(ind[2]^aind[2])*dtz

    txout = txin+dtx
    tyout = tyin+dty
    tzout = tzin+dtz

    tmin = max(txin, tyin, tzin, t0)
    tmax = min(txout, tyout, tzout)

    # Loop over all cells until reaching the out face
    while tmax < tmax_domain:
        octree.get_root(ind, &oct)

        if oct != NULL and tmax > 0:  # no need to check tmin < tmax, as dtx,dty,dtz > 0
            # Hits, so process subtree
            nextDom = proc_subtree(txin, tyin, tzin, txout, tyout, tzout, tmin,
                oct, a, octList, cellList, tList, curDom)
        if tmax == txout:    # Next x
            ind[0] += idx[0]
            txin = txout
            txout += dtx
        elif tmax == tyout:  # Next y
            ind[1] += idx[1]
            tyin = tyout
            tyout += dty
        elif tmax == tzout:  # Next z
            ind[2] += idx[2]
            tzin = tzout
            tzout += dtz

        # Local in/out ts
        tmax = min(txout, tyout, tzout)
        if tList.size() > 0:
            tmin = tList.back()
        else:
            tmin = max(txin, tyin, tzin, t0)

        # When coming from other domain, this is not necessarily true so we want at least one step.
        if nextDom != curDom:
            break

    if tmax >= tmax_domain:  # Return -1 to stop if reached the right edge
        return -1
    else:                    # Return next domain otherwise
        return nextDom

cdef np.uint8_t find_firstNode(
        const np.float64_t tx0, const np.float64_t ty0, const np.float64_t tz0,
        const np.float64_t txM, const np.float64_t tyM, const np.float64_t tzM) nogil:

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
cdef inline np.uint8_t next_node(const np.uint8_t currNode,
        const np.float64_t tx1, const np.float64_t ty1, const np.float64_t tz1) nogil:
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

# Checks if a *cell* is a leaf: either its parent oct has no children, or the cell is not an oct itself
cdef inline bool isLeaf(const Oct *o, const np.uint8_t currNode) nogil:
    return (o.children == NULL) or (o.children[currNode] == NULL)

cdef inline np.uint8_t swap3bits(const np.uint8_t lev) nogil:
    # Swap little-endian to big-endian (C to F ordering)
    return ((lev & 0b100) >> 2) + (lev & 0b010) + ((lev & 0b001) << 2)

# TODO: support negative directions
# TODO: support ray length
@cython.cdivision(True)
@cython.boundscheck(False)
cdef int proc_subtree(
        const np.float64_t tx0, const np.float64_t ty0, const np.float64_t tz0,
        const np.float64_t tx1, const np.float64_t ty1, const np.float64_t tz1,
        const np.float64_t t0,
        const Oct* oct, const int a, vector[np.uint64_t] &octList, vector[np.uint8_t] &cellList, vector[np.float64_t] &tList,
        const int curDom, int level=0) nogil:
    cdef np.uint8_t currNode, nextNode
    cdef np.float64_t txM, tyM, tzM
    cdef np.uint8_t entry_plane, exit_plane
    cdef bool leaf
    cdef int nextDom

    if oct == NULL:
        raise Exception('proc_subtree received a NULL oct. This is likely a bug.')
    nextDom = oct.domain

    if tx1 < t0 or ty1 < t0 or tz1 < t0:
        return oct.domain

    # Note: we cannot break early if the domain is the wrong one, as leaf cells may point to yet another domain
    # so we have to explore the octree down to leaf cells and find the first *leaf* that's not in the current domain

    # Compute midpoints
    txM = (tx0 + tx1) / 2.
    tyM = (ty0 + ty1) / 2.
    tzM = (tz0 + tz1) / 2.

    currNode = find_firstNode(tx0, ty0, tz0, txM, tyM, tzM)
    nextNode = currNode

    while nextNode < 8:
        currNode = nextNode
        leaf = isLeaf(oct, currNode^a)
        # print('%scurrNode=%s %s (%.2f %.2f %.2f) (%.2f %.2f %.2f)' % ('\t'*level, currNode, a, txM, tyM, tzM, tx1, ty1, tz1))
        # Note: there is a bit of code repetition down there (nextNode = ...) but couldn't find a clever way that also efficient
        if leaf:
            # Need to swap bits before storing as octree is C-style in memory and F-style on file
            if curDom == nextDom:
                octList.push_back(oct.domain_ind)
                cellList.push_back(swap3bits(currNode^a))

            if currNode == 0:
                if curDom == nextDom:
                    tList.push_back(min(txM, tyM, tzM))
                nextNode = next_node(currNode, txM, tyM, tzM)
            elif currNode == 1:
                if curDom == nextDom:
                    tList.push_back(min(txM, tyM, tz1))
                nextNode = next_node(currNode, txM, tyM, tz1)
            elif currNode == 2:
                if curDom == nextDom:
                    tList.push_back(min(txM, ty1, tzM))
                nextNode = next_node(currNode, txM, ty1, tzM)
            elif currNode == 3:
                if curDom == nextDom:
                    tList.push_back(min(txM, ty1, tz1))
                nextNode = next_node(currNode, txM, ty1, tz1)
            elif currNode == 4:
                if curDom == nextDom:
                    tList.push_back(min(tx1, tyM, tzM))
                nextNode = next_node(currNode, tx1, tyM, tzM)
            elif currNode == 5:
                if curDom == nextDom:
                    tList.push_back(min(tx1, tyM, tz1))
                nextNode = next_node(currNode, tx1, tyM, tz1)
            elif currNode == 6:
                if curDom == nextDom:
                    tList.push_back(min(tx1, ty1, tzM))
                nextNode = next_node(currNode, tx1, ty1, tzM)
            else:#currNode == 7:
                if curDom == nextDom:
                    tList.push_back(min(tx1, ty1, tz1))
                nextNode = 8

        else:  # Go down the tree
            if currNode == 0:
                nextDom = proc_subtree(tx0, ty0, tz0, txM, tyM, tzM, t0, oct.children[  a], a, octList, cellList, tList, curDom, level+1)
                nextNode = next_node(currNode, txM, tyM, tzM)
            elif currNode == 1:
                nextDom = proc_subtree(tx0, ty0, tzM, txM, tyM, tz1, t0, oct.children[1^a], a, octList, cellList, tList, curDom, level+1)
                nextNode = next_node(currNode, txM, tyM, tz1)
            elif currNode == 2:
                nextDom = proc_subtree(tx0, tyM, tz0, txM, ty1, tzM, t0, oct.children[2^a], a, octList, cellList, tList, curDom, level+1)
                nextNode = next_node(currNode, txM, ty1, tzM)
            elif currNode == 3:
                nextDom = proc_subtree(tx0, tyM, tzM, txM, ty1, tz1, t0, oct.children[3^a], a, octList, cellList, tList, curDom, level+1)
                nextNode = next_node(currNode, txM, ty1, tz1)
            elif currNode == 4:
                nextDom = proc_subtree(txM, ty0, tz0, tx1, tyM, tzM, t0, oct.children[4^a], a, octList, cellList, tList, curDom, level+1)
                nextNode = next_node(currNode, tx1, tyM, tzM)
            elif currNode == 5:
                nextDom = proc_subtree(txM, ty0, tzM, tx1, tyM, tz1, t0, oct.children[5^a], a, octList, cellList, tList, curDom, level+1)
                nextNode = next_node(currNode, tx1, tyM, tz1)
            elif currNode == 6:
                nextDom = proc_subtree(txM, tyM, tz0, tx1, ty1, tzM, t0, oct.children[6^a], a, octList, cellList, tList, curDom, level+1)
                nextNode = next_node(currNode, tx1, ty1, tzM)
            else:#currNode == 7:
                nextDom = proc_subtree(txM, tyM, tzM, tx1, ty1, tz1, t0, oct.children[7  ], a, octList, cellList, tList, curDom, level+1)
                nextNode = 8

    return nextDom

cpdef domain2ind(SparseOctreeContainer octree, SelectorObject selector,
                 np.ndarray[np.int64_t, ndim=2] cell_inds_in, int domain_id=-1):
    cdef StoreIndex visitor
    cdef np.ndarray[np.int64_t, ndim=4] cell_inds

    cell_inds = cell_inds_in.reshape(octree.nocts, 2, 2, 2)

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
