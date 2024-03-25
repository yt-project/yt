cimport numpy as np


cdef (np.float64_t, np.float64_t, np.float64_t) spherical_to_cartesian(np.float64_t r,
                           np.float64_t theta,
                           np.float64_t phi) noexcept nogil


cdef (np.float64_t, np.float64_t, np.float64_t) cartesian_to_spherical(np.float64_t x,
                           np.float64_t y,
                           np.float64_t z) noexcept nogil

cdef class MixedCoordBBox:
    cdef void get_cartesian_bbox(self,
                                np.float64_t pos0,
                                np.float64_t pos1,
                                np.float64_t pos2,
                                np.float64_t dpos0,
                                np.float64_t dpos1,
                                np.float64_t dpos2,
                                np.float64_t xyz_i[3],
                                np.float64_t dxyz_i[3]
                                ) noexcept nogil


cdef class SphericalMixedCoordBBox(MixedCoordBBox):
    cdef void get_cartesian_bbox(
                        self,
                        np.float64_t pos0,
                        np.float64_t pos1,
                        np.float64_t pos2,
                        np.float64_t dpos0,
                        np.float64_t dpos1,
                        np.float64_t dpos2,
                        np.float64_t xyz_i[3],
                        np.float64_t dxyz_i[3]
                        ) noexcept nogil
