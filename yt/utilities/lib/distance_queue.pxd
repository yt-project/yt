"""
A queue for evaluating distances to discrete points




"""


cimport cython
cimport numpy as np

import numpy as np

from libc.stdlib cimport free, malloc
from libc.string cimport memmove

# THESE TWO STRUCTS MUST BE EQUIVALENT

cdef struct ItemList:
    np.int64_t ind
    np.float64_t value

cdef struct NeighborList:
    np.int64_t pn       # Particle number
    np.float64_t r2     # radius**2

cdef int Neighbor_compare(void *on1, void *on2) nogil
cdef np.float64_t r2dist(np.float64_t ppos[3],
                         np.float64_t cpos[3],
                         np.float64_t DW[3],
                         bint periodicity[3],
                         np.float64_t max_dist2)

cdef class PriorityQueue:
    cdef int maxn
    cdef int curn
    cdef ItemList* items
    cdef void item_reset(self)
    cdef int item_insert(self, np.int64_t i, np.float64_t value)

cdef class DistanceQueue(PriorityQueue):
    cdef np.float64_t DW[3]
    cdef bint periodicity[3]
    cdef NeighborList* neighbors # flat array
    cdef void _setup(self, np.float64_t DW[3], bint periodicity[3])
    cdef void neighbor_eval(self, np.int64_t pn, np.float64_t ppos[3],
                            np.float64_t cpos[3])
    cdef void neighbor_reset(self)
