cimport numpy as np

from yt.geometry.oct_container cimport SparseOctreeContainer


cdef class Ray(object):
    cdef np.float64_t origin[3]
    cdef np.float64_t direction[3]
    cdef np.float64_t length

    cpdef np.ndarray[np.float64_t, ndim=1] at(self, np.float64_t t)