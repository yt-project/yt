cimport cython 
import numpy as np
cimport numpy as np
from libc.math cimport sqrt


ctypedef fused Real:
    np.float32_t
    np.float64_t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline Real dot(const Real[3] a, 
                     const Real[3] b) nogil:
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] 


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void cross(const Real[3] a, 
                       const Real[3] b,
                       Real c[3]) nogil:
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void subtract(const Real[3] a, 
                          const Real[3] b,
                          Real c[3]) nogil:
    c[0] = a[0] - b[0]
    c[1] = a[1] - b[1]
    c[2] = a[2] - b[2]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void fma(const Real f,
                     const Real[3] a, 
                     const Real[3] b,
                     Real[3] c) nogil:
    c[0] = f * a[0] + b[0]
    c[1] = f * a[1] + b[1]
    c[2] = f * a[2] + b[2]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline Real L2_norm(const Real[3] a) nogil:
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
