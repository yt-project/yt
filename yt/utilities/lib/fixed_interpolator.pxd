"""
Fixed interpolator includes



"""


cimport numpy as np


cdef extern from "fixed_interpolator.hpp":
    np.float64_t fast_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                  np.float64_t *data) noexcept nogil
    np.float64_t offset_interpolate(int ds[3], np.float64_t dp[3],
                                    np.float64_t *data) noexcept nogil
    np.float64_t trilinear_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                       np.float64_t *data) noexcept nogil
    void eval_gradient(int ds[3], np.float64_t dp[3], np.float64_t *data,
                       np.float64_t grad[3]) noexcept nogil
    void offset_fill(int *ds, np.float64_t *data, np.float64_t *gridval) noexcept nogil
    void vertex_interp(np.float64_t v1, np.float64_t v2, np.float64_t isovalue,
                       np.float64_t vl[3], np.float64_t dds[3],
                       np.float64_t x, np.float64_t y, np.float64_t z,
                       int vind1, int vind2) noexcept nogil
