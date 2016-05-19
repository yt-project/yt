"""
Fixed interpolator includes



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np


cdef extern from "fixed_interpolator.h":
    np.float64_t offset_interpolate(int ds[3], np.float64_t dp[3],
                                    np.float64_t *data) nogil
    np.float64_t trilinear_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                       np.float64_t *data) nogil
    void eval_gradient(int ds[3], np.float64_t dp[3], np.float64_t *data,
                       np.float64_t grad[3]) nogil
    void offset_fill(int *ds, np.float64_t *data, np.float64_t *gridval) nogil
