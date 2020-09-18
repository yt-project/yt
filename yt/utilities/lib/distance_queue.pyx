
# distutils: libraries = STD_LIBS
"""
Distance queue implementation




"""

cimport numpy as np

import numpy as np

cimport cython


cdef int Neighbor_compare(void *on1, void *on2) nogil:
    cdef NeighborList *n1
    cdef NeighborList *n2
    n1 = <NeighborList *> on1
    n2 = <NeighborList *> on2
    # Note that we set this up so that "greatest" evaluates to the *end* of the
    # list, so we can do standard radius comparisons.
    if n1.r2 < n2.r2:
        return -1
    elif n1.r2 == n2.r2:
        return 0
    else:
        return 1

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef np.float64_t r2dist(np.float64_t ppos[3],
                         np.float64_t cpos[3],
                         np.float64_t DW[3],
                         bint periodicity[3],
                         np.float64_t max_dist2):
    cdef int i
    cdef np.float64_t r2, DR
    r2 = 0.0
    for i in range(3):
        DR = (ppos[i] - cpos[i])
        if not periodicity[i]:
            pass
        elif (DR > DW[i]/2.0):
            DR -= DW[i]
        elif (DR < -DW[i]/2.0):
            DR += DW[i]
        r2 += DR * DR
        if max_dist2 >= 0.0 and r2 > max_dist2:
            return -1.0
    return r2

cdef class PriorityQueue:
    """This class acts as a "priority-queue."  It was extracted from the
    DistanceQueue object so that it could serve as a general-purpose method for
    storing the N-most "valuable" objects.  It's relatively simple, in that it
    provides storage for a single int64 (which is usually an 'index' into an
    external array or list) and a single float64 "value" associated with them.
    You can insert new objects, and then if it's already at maxn objects in it,
    it'll bump the least valuable one off the end.  Of particular note is that
    *lower* values are considered to be *more valuable* than higher ones; this
    is because our typical use case is to store radii.
    """
    def __cinit__(self, int maxn):
        cdef int i
        self.maxn = maxn
        self.curn = 0
        self.items = <ItemList *> malloc(
            sizeof(ItemList) * self.maxn)

    cdef void item_reset(self):
        for i in range(self.maxn):
            self.items[i].value = 1e300
            self.items[i].ind = -1
        self.curn = 0

    cdef int item_insert(self, np.int64_t ind, np.float64_t value):
        cdef int i, di
        if self.curn == 0:
            self.items[0].value = value
            self.items[0].ind = ind
            self.curn += 1
            return 0
        # Now insert in a sorted way
        di = -1
        for i in range(self.curn - 1, -1, -1):
            # We are checking if i is less than us, to see if we should insert
            # to the right (i.e., i+1).
            if self.items[i].value < value:
                di = i
                break
        # The outermost one is already too small.
        if di == self.maxn - 1:
            return -1
        if (self.maxn - (di + 2)) > 0:
            memmove(<void *> (self.items + di + 2),
                    <void *> (self.items + di + 1),
                    sizeof(ItemList) * (self.maxn - (di + 2)))
        self.items[di + 1].value = value
        self.items[di + 1].ind = ind
        if self.curn < self.maxn:
            self.curn += 1
        return di + 1

cdef class DistanceQueue:
    """This is a distance queue object, designed to incrementally evaluate N
    positions against a reference point.  It is initialized with the maximum
    number that are to be retained (i.e., maxn-nearest neighbors)."""
    def __cinit__(self, int maxn):
        if sizeof(ItemList) != sizeof(NeighborList):
            # This should almost never, ever happen unless we do something very
            # wrong, and must be broken at compile time.
            raise RuntimeError
        self.neighbors = <NeighborList *> self.items
        self.neighbor_reset()
        for i in range(3):
            self.DW[i] = 0
            self.periodicity[i] = False

    cdef void _setup(self, np.float64_t DW[3], bint periodicity[3]):
        for i in range(3):
            self.DW[i] = DW[i]
            self.periodicity[i] = periodicity[i]

    def setup(self, np.float64_t[:] DW, periodicity):
        for i in range(3):
            self.DW[i] = DW[i]
            self.periodicity[i] = periodicity[i]

    def __dealloc__(self):
        free(self.neighbors)

    cdef void neighbor_eval(self, np.int64_t pn, np.float64_t ppos[3],
                            np.float64_t cpos[3]):
        # Here's a python+numpy simulator of this:
        # http://paste.yt-project.org/show/5445/
        cdef np.float64_t r2, r2_trunc
        if self.curn == self.maxn:
            # Truncate calculation if it's bigger than this in any dimension
            r2_trunc = self.neighbors[self.curn - 1].r2
        else:
            # Don't truncate our calculation
            r2_trunc = -1
        r2 = r2dist(ppos, cpos, self.DW, self.periodicity, r2_trunc)
        if r2 == -1:
            return
        self.item_insert(pn, r2)

    cdef void neighbor_reset(self):
        self.item_reset()

    def find_nearest(self, np.float64_t[:] center, np.float64_t[:,:] points):
        """This function accepts a center and a set of [N,3] points, and it
        will return the indices into the points array of the maxn closest
        neighbors."""
        cdef int i, j
        cdef np.float64_t ppos[3]
        cdef np.float64_t cpos[3]
        self.neighbor_reset()
        for i in range(3):
            cpos[i] = center[i]
        for j in range(points.shape[0]):
            for i in range(3):
                ppos[i] = points[j,i]
            self.neighbor_eval(j, ppos, cpos)
        rv = np.empty(self.curn, dtype="int64")
        for i in range(self.curn):
            rv[i] = self.neighbors[i].pn
        return rv
