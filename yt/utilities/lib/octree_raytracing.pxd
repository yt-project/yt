cimport numpy as np
from libcpp.vector cimport vector
from libc.stdlib cimport bsearch, qsort, realloc, malloc, free

from yt.geometry.oct_container cimport SparseOctreeContainer


cdef struct OctreeDescriptor:
    # SparseOctreeContainer octree
    np.float64_t* oct_LE
    np.float64_t* oct_ddd
    int nocts
    int domain_id
    int nfields
    np.float64_t* data

cdef class OctreesDescriptor:
    cdef OctreeDescriptor *oct_descs
    cdef dict octrees
    cdef dict dom2ind

    cdef OctreeDescriptor *get_descriptor(self, int domain_id) nogil

cdef int DomainDecomposition_find_domain(np.float64_t bscale, int bit_length, int ncpu, np.uint64_t[:] keys, np.float64_t[:] pos) nogil

cdef int ray_step(SparseOctreeContainer octree,
                  np.float64_t* origin, np.float64_t* direction,
                  vector[np.uint64_t] &octList, vector[np.uint8_t] &cellList, vector[np.float64_t] &tList,
                  np.float64_t t0, const int curDom) nogil
