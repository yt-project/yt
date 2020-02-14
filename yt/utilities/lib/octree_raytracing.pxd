cimport numpy as np
from libcpp.vector cimport vector

from yt.geometry.oct_container cimport SparseOctreeContainer

cdef int DomainDecomposition_find_domain(np.float64_t bscale, int bit_length, int ncpu, np.uint64_t[:] keys, np.float64_t[:] pos) nogil

cdef int ray_step(SparseOctreeContainer octree,
                  np.float64_t* origin, np.float64_t* direction,
                  vector[np.uint64_t] &octList, vector[np.uint8_t] &cellList, vector[np.float64_t] &tList,
                  np.float64_t t0, const int curDom) nogil
