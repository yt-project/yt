cimport cython 
cimport cython.floating
from libc.math cimport sqrt


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline cython.floating dot(const cython.floating[3] a, 
                                const cython.floating[3] b) nogil:
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] 


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void cross(const cython.floating[3] a,
                       const cython.floating[3] b,
                       cython.floating c[3]) nogil:
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void subtract(const cython.floating[3] a, 
                          const cython.floating[3] b,
                          cython.floating c[3]) nogil:
    c[0] = a[0] - b[0]
    c[1] = a[1] - b[1]
    c[2] = a[2] - b[2]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline cython.floating distance(const cython.floating[3] a,
                                     const cython.floating[3] b) nogil:
    return sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2 +(a[2] - b[2])**2)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void fma(const cython.floating f,
                     const cython.floating[3] a, 
                     const cython.floating[3] b,
                     cython.floating[3] c) nogil:
    c[0] = f * a[0] + b[0]
    c[1] = f * a[1] + b[1]
    c[2] = f * a[2] + b[2]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline cython.floating L2_norm(const cython.floating[3] a) nogil:
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
