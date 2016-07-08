"""
A queue for evaluating distances to discrete points




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport cython
cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free
from libc.string cimport memmove

cdef struct NeighborList:
    np.int64_t pn       # Particle number
    np.float64_t r2     # radius**2

cdef int Neighbor_compare(void *on1, void *on2) nogil
cdef np.float64_t r2dist(np.float64_t ppos[3],
                         np.float64_t cpos[3],
                         np.float64_t DW[3],
                         bint periodicity[3],
                         np.float64_t max_dist2)

cdef class DistanceQueue:
    cdef int maxn
    cdef int curn
    cdef np.float64_t DW[3]
    cdef bint periodicity[3]
    cdef NeighborList* neighbors # flat array
    cdef void _setup(self, np.float64_t DW[3], bint periodicity[3])
    cdef void neighbor_eval(self, np.int64_t pn, np.float64_t ppos[3],
                            np.float64_t cpos[3])
    cdef void neighbor_reset(self)
