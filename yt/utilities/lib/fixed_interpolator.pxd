# cython: cdivision=True
"""
Fixed interpolator includes



"""


cimport numpy as np
from libc.math cimport fabs, sqrt


cdef inline int VINDEX(int A, int B, int C, int ds[3], int ci[3]) nogil:
    return ((((A)+ci[0])*(ds[1]+1)+((B)+ci[1]))*(ds[2]+1)+ci[2]+(C))

cdef inline int OINDEX(int A, int B, int C, int ds[3]) nogil:
    return (A)*(ds[1]+1)*(ds[2]+1)+(B)*ds[2]+(B)+(C)

cdef inline np.float64_t fast_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                   np.float64_t *data) nogil:
    cdef int i
    cdef np.float64_t dv
    cdef np.float64_t dm[3]
    for i in range(3):
        dm[i] = (1.0 - dp[i])
    dv  = 0.0;
    dv += data[VINDEX(0,0,0, ds, ci)] * (dm[0]*dm[1]*dm[2])
    dv += data[VINDEX(0,0,1, ds, ci)] * (dm[0]*dm[1]*dp[2])
    dv += data[VINDEX(0,1,0, ds, ci)] * (dm[0]*dp[1]*dm[2])
    dv += data[VINDEX(0,1,1, ds, ci)] * (dm[0]*dp[1]*dp[2])
    dv += data[VINDEX(1,0,0, ds, ci)] * (dp[0]*dm[1]*dm[2])
    dv += data[VINDEX(1,0,1, ds, ci)] * (dp[0]*dm[1]*dp[2])
    dv += data[VINDEX(1,1,0, ds, ci)] * (dp[0]*dp[1]*dm[2])
    dv += data[VINDEX(1,1,1, ds, ci)] * (dp[0]*dp[1]*dp[2])
    return dv

cdef inline np.float64_t offset_interpolate(int ds[3], np.float64_t dp[3],
                                     np.float64_t *data) nogil:
    cdef int i
    cdef np.float64_t dv
    cdef np.float64_t vz[4]

    dv = 1.0 - dp[2]
    vz[0] = dv*data[OINDEX(0,0,0,ds)] + dp[2]*data[OINDEX(0,0,1,ds)]
    vz[1] = dv*data[OINDEX(0,1,0,ds)] + dp[2]*data[OINDEX(0,1,1,ds)]
    vz[2] = dv*data[OINDEX(1,0,0,ds)] + dp[2]*data[OINDEX(1,0,1,ds)]
    vz[3] = dv*data[OINDEX(1,1,0,ds)] + dp[2]*data[OINDEX(1,1,1,ds)]

    dv = 1.0 - dp[1]
    vz[0] = dv*vz[0] + dp[1]*vz[1]
    vz[1] = dv*vz[2] + dp[1]*vz[3]

    dv = 1.0 - dp[0]
    vz[0] = dv*vz[0] + dp[0]*vz[1]

    return vz[0]

cdef inline np.float64_t trilinear_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                       np.float64_t *data) nogil:
    cdef int i
    cdef np.float64_t dm[3]
    cdef np.float64_t vz[4]
    # dp is the distance to the plane.  dm is val, dp = 1-val
    for i in range(3):
        dm[i] = 1.0 - dp[i]

    # First interpolate in z
    vz[0] = dm[2]*data[VINDEX(0,0,0,ds,ci)] + dp[2]*data[VINDEX(0,0,1,ds,ci)]
    vz[1] = dm[2]*data[VINDEX(0,1,0,ds,ci)] + dp[2]*data[VINDEX(0,1,1,ds,ci)]
    vz[2] = dm[2]*data[VINDEX(1,0,0,ds,ci)] + dp[2]*data[VINDEX(1,0,1,ds,ci)]
    vz[3] = dm[2]*data[VINDEX(1,1,0,ds,ci)] + dp[2]*data[VINDEX(1,1,1,ds,ci)]

    # Then in y
    vz[0] = dm[1]*vz[0] + dp[1]*vz[1]
    vz[1] = dm[1]*vz[2] + dp[1]*vz[3]

    # Then in x
    vz[0] = dm[0]*vz[0] + dp[0]*vz[1];
    return vz[0]


cdef inline void eval_gradient(int ds[3], np.float64_t dp[3], np.float64_t *data,
                               np.float64_t grad[3]) nogil:
    # We just take some small value

    cdef int i
    cdef np.float64_t denom, plus, minus, backup, normval

    normval = 0.0;
    for i in range(3):
        backup = dp[i]
        grad[i] = 0.0
        if (dp[i] >= 0.95):
            plus = dp[i]
            minus = dp[i] - 0.05
        elif dp[i] <= 0.05:
            plus = dp[i] + 0.05
            minus = 0.0
        else:
            plus = dp[i] + 0.05
            minus = dp[i] - 0.05
        denom = plus - minus
        dp[i] = plus
        grad[i] += offset_interpolate(ds, dp, data) / denom
        dp[i] = minus
        grad[i] -= offset_interpolate(ds, dp, data) / denom
        dp[i] = backup
        normval += grad[i]*grad[i]
    if normval != 0.0:
        normval = sqrt(normval)
        for i in range(3):
            grad[i] /= -normval
    else:
        grad[0] = grad[1] = grad[2] = 0.0

cdef inline void offset_fill(int *ds, np.float64_t *data, np.float64_t *gridval) nogil:
    gridval[0] = data[OINDEX(0,0,0,ds)]
    gridval[1] = data[OINDEX(1,0,0,ds)]
    gridval[2] = data[OINDEX(1,1,0,ds)]
    gridval[3] = data[OINDEX(0,1,0,ds)]
    gridval[4] = data[OINDEX(0,0,1,ds)]
    gridval[5] = data[OINDEX(1,0,1,ds)]
    gridval[6] = data[OINDEX(1,1,1,ds)]
    gridval[7] = data[OINDEX(0,1,1,ds)]


cdef inline void vertex_interp(np.float64_t v1, np.float64_t v2, np.float64_t isovalue,
                       np.float64_t vl[3], np.float64_t dds[3],
                       np.float64_t x, np.float64_t y, np.float64_t z,
                               int vind1, int vind2) nogil:
    cdef int i
    cdef np.float64_t mu = ((isovalue - v1) / (v2 - v1))
    cdef np.float64_t* cverts[8]
    cverts[0] = [0, 0, 0]
    cverts[1] = [1, 0, 0]
    cverts[2] = [1, 1, 0]
    cverts[3] = [0, 1, 0]
    cverts[4] = [0, 0, 1]
    cverts[5] = [1, 0, 1]
    cverts[6] = [1, 1, 1]
    cverts[7] = [0, 1, 1]

    if (fabs(1.0 - isovalue/v1) < 0.000001):
        mu = 0.0

    if (fabs(1.0 - isovalue/v2) < 0.000001):
        mu = 1.0
    if (fabs(v1/v2) < 0.000001):
        mu = 0.0;

    vl[0] = x
    vl[1] = y
    vl[2] = z

    for i in range(3):
        vl[i] += (dds[i] * cverts[vind1][i]
               + dds[i] * mu*(cverts[vind2][i] - cverts[vind1][i]))
