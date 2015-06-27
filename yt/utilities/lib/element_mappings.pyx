"""
This file contains coordinate mappings between physical coordinates and those
defined on unit elements, as well as doing the corresponding intracell 
interpolation on finite element data.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np
from numpy cimport ndarray
cimport cython
import numpy as np
from libc.math cimport abs


cdef class ElementSampler:

    def map_real_to_unit(self,
                         np.ndarray physical_coord, 
                         np.ndarray vertices):
        raise NotImplementedError

    def sample_at_unit_point(self,
                             np.ndarray coord, 
                             np.ndarray vals):
        raise NotImplementedError

    def sample_at_real_point(self,
                             np.ndarray coord, 
                             np.ndarray vertices, 
                             np.ndarray vals):
        mapped_coord = self.map_real_to_unit(coord, vertices)
        return self.sample_at_unit_point(mapped_coord, vals)
    

cdef class P1Sampler2D(ElementSampler):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def map_real_to_unit(self, 
                         np.ndarray[np.float64_t, ndim=1] physical_coord, 
                         np.ndarray[np.float64_t, ndim=2] vertices):
    
        x = physical_coord[0]
        y = physical_coord[1]

        x1 = vertices[0, 0]
        y1 = vertices[0, 1]

        x2 = vertices[1, 0]
        y2 = vertices[1, 1]

        x3 = vertices[2, 0]
        y3 = vertices[2, 1]
    
        A = np.array([[1, x, y], [1, x1, y1], [1, x3, y3]])
        B = np.array([[1, x2, y2], [1, x1, y1], [1, x3, y3]])
        u = np.linalg.det(A) / np.linalg.det(B)

        C = np.array([[1, x, y], [1, x1, y1], [1, x2, y2]])
        D = np.array([[1, x3, y3], [1, x1, y1], [1, x2, y2]])
        v = np.linalg.det(C) / np.linalg.det(D)
    
        return np.array([u, v])

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def sample_at_unit_point(self,
                             np.ndarray coord, 
                             np.ndarray vals):
        return vals[0]*(1 - coord[0] - coord[1]) + \
            vals[1]*coord[0] + vals[2]*coord[1]

cdef class P1Sampler3D(ElementSampler):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def map_real_to_unit(self, 
                         np.ndarray[np.float64_t, ndim=1] physical_coord, 
                         np.ndarray[np.float64_t, ndim=2] vertices):
    
        b = np.array([physical_coord[0],
                      physical_coord[1], 
                      physical_coord[2], 
                      1.0], dtype=np.float64)

        A = np.empty((4, 4), dtype=np.float64)
        cdef int i, j
        for i in range(3):
            for j in range(4):
                A[i][j] = vertices[j][i]
        for j in range(4):
            A[3][j] = 1.0

        c = np.linalg.solve(A, b)
    
        return c

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def sample_at_unit_point(self,
                             np.ndarray[np.float64_t, ndim=1] coord, 
                             np.ndarray[np.float64_t, ndim=1] vals):
        cdef np.float64_t value = 0.0
        cdef np.int64_t i
        for i in range(4):
            value += vals[i]*coord[i]
        return value

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef inline void Q1Function2D(double[:] fx,
                               double[:] x, 
                               double[:, :] vertices, 
                               double[:] phys_x) nogil:
    
    cdef int i
    for i in range(2):
        fx[i] = vertices[0][i]*(1-x[0])*(1-x[1]) \
              + vertices[1][i]*(1+x[0])*(1-x[1]) \
              + vertices[2][i]*(1-x[0])*(1+x[1]) \
              + vertices[3][i]*(1+x[0])*(1+x[1]) - 4.0*phys_x[i]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef inline void Q1Jacobian2D(double[:, :] A,
                               double[:] x, 
                               double[:, :] v, 
                               double[:] phys_x) nogil:
    
    cdef int i
    for i in range(2):
        A[i][0] = -(1-x[1])*v[0][i] + (1-x[1])*v[1][i] - \
                   (1+x[1])*v[2][i] + (1+x[1])*v[3][i]
        A[i][1] = -(1-x[0])*v[0][i] - (1+x[0])*v[1][i] + \
                   (1-x[0])*v[2][i] + (1+x[0])*v[3][i]


cdef class Q1Sampler2D(ElementSampler):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def map_real_to_unit(self, 
                         np.ndarray[np.float64_t, ndim=1] physical_x,
                         np.ndarray[np.float64_t, ndim=2] vertices):
        x = np.zeros(2, dtype=np.float64)
        cdef int iterations = 0
        cdef np.float64_t tolerance = 1.0e-9
        fx = np.empty(2, dtype=np.float64)
        A = np.empty((2, 2), dtype=np.float64)
        Ainv = np.empty((2, 2), dtype=np.float64)
        Q1Function2D(fx, x, vertices, physical_x)
        cdef np.float64_t err = np.max(abs(fx))
        while (err > tolerance and iterations < 100):
            Q1Jacobian2D(A, x, vertices, physical_x)
            Ainv = np.linalg.inv(A)
            x = x - np.dot(Ainv, fx)
            Q1Function2D(fx, x, vertices, physical_x)
            err = np.max(abs(fx))
            iterations += 1
        return x

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def sample_at_unit_point(self, double[:] coord, 
                             double[:] vals):
        cdef double x = vals[0]*(1.0 - coord[0])*(1.0 - coord[1]) + \
                        vals[1]*(1.0 + coord[0])*(1.0 - coord[1]) + \
                        vals[2]*(1.0 - coord[0])*(1.0 + coord[1]) + \
                        vals[3]*(1.0 + coord[0])*(1.0 + coord[1])
        return 0.25*x

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef inline void Q1Function3D(double[:] fx,
                               double[:] x, 
                               double[:, :] vertices, 
                               double[:] phys_x) nogil:
    
    cdef int i
    for i in range(3):
        fx[i] = vertices[0][i]*(1-x[0])*(1-x[1])*(1-x[2]) \
              + vertices[1][i]*(1+x[0])*(1-x[1])*(1-x[2]) \
              + vertices[2][i]*(1-x[0])*(1+x[1])*(1-x[2]) \
              + vertices[3][i]*(1+x[0])*(1+x[1])*(1-x[2]) \
              + vertices[4][i]*(1-x[0])*(1-x[1])*(1+x[2]) \
              + vertices[5][i]*(1+x[0])*(1-x[1])*(1+x[2]) \
              + vertices[6][i]*(1-x[0])*(1+x[1])*(1+x[2]) \
              + vertices[7][i]*(1+x[0])*(1+x[1])*(1+x[2]) \
              - 8.0*phys_x[i]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef inline void Q1Jacobian3D(double[:, :] A,
                               double[:] x, 
                               double[:, :] v, 
                               double[:] phys_x) nogil:
    
    cdef int i
    for i in range(3):
        A[i][0] = -(1-x[1])*(1-x[2])*v[0][i] + (1-x[1])*(1-x[2])*v[1][i] - \
                   (1+x[1])*(1-x[2])*v[2][i] + (1+x[1])*(1-x[2])*v[3][i] - \
                   (1-x[1])*(1+x[2])*v[4][i] + (1-x[1])*(1+x[2])*v[5][i] - \
                   (1+x[1])*(1+x[2])*v[6][i] + (1+x[1])*(1+x[2])*v[7][i]
        A[i][1] = -(1-x[0])*(1-x[2])*v[0][i] - (1+x[0])*(1-x[2])*v[1][i] + \
                   (1-x[0])*(1-x[2])*v[2][i] + (1+x[0])*(1-x[2])*v[3][i] - \
                   (1-x[0])*(1+x[2])*v[4][i] - (1+x[0])*(1+x[2])*v[5][i] + \
                   (1-x[0])*(1+x[2])*v[6][i] + (1+x[0])*(1+x[2])*v[7][i]
        A[i][2] = -(1-x[0])*(1-x[1])*v[0][i] - (1+x[0])*(1-x[1])*v[1][i] - \
                   (1-x[0])*(1+x[1])*v[2][i] - (1+x[0])*(1+x[1])*v[3][i] + \
                   (1-x[0])*(1-x[1])*v[4][i] + (1+x[0])*(1-x[1])*v[5][i] + \
                   (1-x[0])*(1+x[1])*v[6][i] + (1+x[0])*(1+x[1])*v[7][i]


cdef class Q1Sampler3D(ElementSampler):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def map_real_to_unit(self, 
                         np.ndarray[np.float64_t, ndim=1] physical_x,
                         np.ndarray[np.float64_t, ndim=2] vertices):
        x = np.zeros(3, dtype=np.float64)
        cdef int iterations = 0
        cdef np.float64_t tolerance = 1.0e-9
        fx = np.empty(3, dtype=np.float64)
        A = np.empty((3, 3), dtype=np.float64)
        Ainv = np.empty((3, 3), dtype=np.float64)
        Q1Function3D(fx, x, vertices, physical_x)
        cdef np.float64_t err = np.max(abs(fx))
        while (err > tolerance and iterations < 100):
            Q1Jacobian3D(A, x, vertices, physical_x)
            Ainv = np.linalg.inv(A)
            x = x - np.dot(Ainv, fx)
            Q1Function3D(fx, x, vertices, physical_x)
            err = np.max(abs(fx))
            iterations += 1
        return x

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def sample_at_unit_point(self, double[:] coord, double[:] vals):
        cdef double x = vals[0]*(1.0 - coord[0])*(1.0 - coord[1])*(1.0 - coord[2]) + \
                        vals[1]*(1.0 + coord[0])*(1.0 - coord[1])*(1.0 - coord[2]) + \
                        vals[2]*(1.0 - coord[0])*(1.0 + coord[1])*(1.0 - coord[2]) + \
                        vals[3]*(1.0 + coord[0])*(1.0 + coord[1])*(1.0 - coord[2]) + \
                        vals[4]*(1.0 - coord[0])*(1.0 - coord[1])*(1.0 + coord[2]) + \
                        vals[5]*(1.0 + coord[0])*(1.0 - coord[1])*(1.0 + coord[2]) + \
                        vals[6]*(1.0 - coord[0])*(1.0 + coord[1])*(1.0 + coord[2]) + \
                        vals[7]*(1.0 + coord[0])*(1.0 + coord[1])*(1.0 + coord[2])
        return 0.125*x
