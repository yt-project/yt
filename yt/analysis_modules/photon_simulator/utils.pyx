import numpy as np
cimport numpy as np
cimport cython

cdef extern from "platform_dep.h":
    double erf(double x)
    
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def broaden_lines(np.ndarray[np.float64_t, ndim=1] E0,
                  np.ndarray[np.float64_t, ndim=1] sigma,
                  np.ndarray[np.float64_t, ndim=1] amp,
                  np.ndarray[np.float64_t, ndim=1] ebins):

    cdef int i, j, n, m
    cdef double x, isigma
    cdef np.ndarray[np.float64_t, ndim=1] cdf, vec

    n = E0.shape[0]
    m = ebins.shape[0]
    cdf = np.zeros(m)
    vec = np.zeros(m-1)
    
    for i in range(n):
        isigma = 1.0/sigma[i]
        for j in range(m):
            x = (ebins[j]-E0[i])*isigma
            cdf[j] = 0.5*(1+erf(x))
        for j in range(m-1):
            vec[j] = vec[j] + (cdf[j+1] - cdf[j])*amp[i]
    return vec
