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

cdef inline np.float64_t offset_interpolate(np.float64_t dp[3],
                np.float64_t[:,:,:] data, int i, int j, int k) nogil:
    cdef np.float64_t dv, vz[4]
    dv = 1.0 - dp[2]
    vz[0] = dv*data[i+0,j+0,k+0] + dp[2]*data[i+0,j+0,k+1]
    vz[1] = dv*data[i+0,j+1,k+0] + dp[2]*data[i+0,j+1,k+1]
    vz[2] = dv*data[i+1,j+0,k+0] + dp[2]*data[i+1,j+0,k+1]
    vz[3] = dv*data[i+1,j+1,k+0] + dp[2]*data[i+1,j+1,k+1]

    dv = 1.0 - dp[1]
    vz[0] = dv*vz[0] + dp[1]*vz[1]
    vz[1] = dv*vz[2] + dp[1]*vz[3]

    dv = 1.0 - dp[0]
    vz[0] = dv*vz[0] + dp[0]*vz[1]

    return vz[0]

cdef inline np.float64_t trilinear_interpolate(int ds[3], int ci[3],
                np.float64_t dp[3], np.float64_t[:,:,:] data) nogil:
    return 0

cdef inline void eval_gradient(int ds[3], np.float64_t dp[3],
                np.float64_t[:,:,] data, np.float64_t grad[3]) nogil:
    return

cdef inline void offset_fill(np.float64_t[:,:,:] data,
                             np.float64_t[:,:,:] gridval,
                             int i, int j, int k) nogil:
    cdef int oi, oj, ok
    for oi in range(2):
        for oj in range(2):
            for ok in range(2):
                data[i+oi, j+oj, k+ok] = gridval[i+oi, j+oj, k+ok]
