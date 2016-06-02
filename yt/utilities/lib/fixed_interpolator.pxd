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

from libc.math cimport sqrt, fabs
cimport cython
cimport numpy as np

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
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

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void eval_gradient(np.float64_t dp[3],
                np.float64_t[:,:,:] data, np.float64_t grad[3],
                int i, int j, int k) nogil:

    # We just take some small value

    cdef int m
    cdef np.float64_t denom, plus, minus, backup, normval
    
    normval = 0.0
    for m in range(3):
        backup = dp[m]
        grad[m] = 0.0
        if (dp[m] >= 0.95):
            plus = dp[m]
            minus = dp[m] - 0.05
        elif (dp[i] <= 0.05):
            plus = dp[m] + 0.05
            minus = 0.0
        else:
            plus = dp[m] + 0.05
            minus = dp[m] - 0.05
        denom = plus - minus
        dp[m] = plus
        grad[m] += offset_interpolate(dp, data, i, j, k) / denom
        dp[m] = minus
        grad[m] -= offset_interpolate(dp, data, i, j, k) / denom
        dp[m] = backup
        normval += grad[m] * grad[m]
    if (normval != 0.0):
        normval = sqrt(normval)
        for m in range(3):
            grad[m] /= -normval
    else:
      grad[0] = grad[1] = grad[2] = 0.0

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void offset_fill(np.float64_t[:,:,:] inval,
                             np.float64_t[:] outval,
                             int i, int j, int k) nogil:
    outval[0] = inval[i+0,j+0,k+0]
    outval[1] = inval[i+1,j+0,k+0]
    outval[2] = inval[i+1,j+1,k+0]
    outval[3] = inval[i+0,j+1,k+0]
    outval[4] = inval[i+0,j+0,k+1]
    outval[5] = inval[i+1,j+0,k+1]
    outval[6] = inval[i+1,j+1,k+1]
    outval[7] = inval[i+0,j+1,k+1]

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void vertex_interp(np.float64_t v1, np.float64_t v2,
            np.float64_t isovalue, np.float64_t vl[3], np.float64_t dds[3],
            np.float64_t x, np.float64_t y, np.float64_t z,
            int vind1, int vind2):
    cdef int i
    cdef np.float64_t* cverts = \
        [0,0,0, 1,0,0, 1,1,0, 0,1,0,
         0,0,1, 1,0,1, 1,1,1, 0,1,1]

    cdef np.float64_t mu = ((isovalue - v1) / (v2 - v1))

    if (fabs(1.0 - isovalue/v1) < 0.000001):
        mu = 0.0
    if (fabs(1.0 - isovalue/v2) < 0.000001):
        mu = 1.0
    if (fabs(v1/v2) < 0.000001):
        mu = 0.0

    vl[0] = x
    vl[1] = y
    vl[2] = z
    for i in range(3):
        vl[i] += dds[i] * cverts[vind1*3+i] \
               + dds[i] * mu*(cverts[vind2*3+i] - cverts[vind1*3+i])
