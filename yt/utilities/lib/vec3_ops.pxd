cimport cython 
import numpy as np
cimport numpy as np
from libc.math cimport sqrt

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline np.float64_t dot(const np.float64_t a[3], 
                             const np.float64_t b[3]) nogil:
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] 


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void cross(const np.float64_t a[3], 
                       const np.float64_t b[3],
                       np.float64_t c[3]) nogil:
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void subtract(const np.float64_t a[3], 
                          const np.float64_t b[3],
                          np.float64_t c[3]) nogil:
    c[0] = a[0] - b[0]
    c[1] = a[1] - b[1]
    c[2] = a[2] - b[2]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void fma(const np.float64_t f,
                     const np.float64_t a[3], 
                     const np.float64_t b[3],
                     np.float64_t c[3]) nogil:
    c[0] = f * a[0] + b[0]
    c[1] = f * a[1] + b[1]
    c[2] = f * a[2] + b[2]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline np.float64_t L2_norm(const np.float64_t a[3]) nogil:
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
