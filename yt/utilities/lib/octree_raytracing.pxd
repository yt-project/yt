cimport numpy as np

from yt.geometry.oct_container cimport SparseOctreeContainer


cdef class Ray(object):
    cdef np.float64_t origin[3]
    cdef np.float64_t direction[3]
    cdef np.float64_t length

    cpdef np.ndarray[np.float64_t, ndim=1] at(self, np.float64_t t)

    cpdef void trilinear(self, const np.float64_t tmin, const np.float64_t tmax,
                         const np.float64_t[:, :, :] data_in, 
                         np.float64_t[:] data_out,
                         const np.float64_t[:] DLE,
                         const np.float64_t Deltax,
                         const int npt)