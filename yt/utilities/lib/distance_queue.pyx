"""
Distance queue implementation




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
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

cdef class DistanceQueue:
    def __cinit__(self, int maxn):
        cdef int i
        self.maxn = maxn
        self.curn = 0
        self.neighbors = <NeighborList *> malloc(
            sizeof(NeighborList) * self.maxn)
        self.neighbor_reset()
        for i in range(3):
            self.DW[i] = 0
            self.periodicit[i] = False

    cdef void setup(self, np.float64_t DW[3], bint periodicity[3]):
        for i in range(3):
            self.DW[i] = DW[i]
            self.periodicity[i] = periodicity[i]

    def __dealloc__(self):
        free(self.neighbors)

    cdef void neighbor_eval(self, np.int64_t pn, np.float64_t ppos[3],
                            np.float64_t cpos[3]):
        # Here's a python+numpy simulator of this:
        # http://paste.yt-project.org/show/5445/
        cdef int i, di
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
        if self.curn == 0:
            self.neighbors[0].r2 = r2
            self.neighbors[0].pn = pn
            self.curn += 1
            return
        # Now insert in a sorted way
        di = -1
        for i in range(self.curn - 1, -1, -1):
            # We are checking if i is less than us, to see if we should insert
            # to the right (i.e., i+1).
            if self.neighbors[i].r2 < r2:
                di = i
                break
        # The outermost one is already too small.
        if di == self.maxn - 1:
            return
        if (self.maxn - (di + 2)) > 0:
            memmove(<void *> (self.neighbors + di + 2),
                    <void *> (self.neighbors + di + 1),
                    sizeof(NeighborList) * (self.maxn - (di + 2)))
        self.neighbors[di + 1].r2 = r2
        self.neighbors[di + 1].pn = pn
        if self.curn < self.maxn:
            self.curn += 1

    cdef void neighbor_reset(self):
        for i in range(self.maxn):
            self.neighbors[i].r2 = 1e300
            self.neighbors[i].pn = -1
        self.curn = 0
