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
    @cython.initializedcheck(False)
    def map_real_to_unit(self, 
                         np.ndarray[np.float64_t, ndim=1] physical_x, 
                         np.ndarray[np.float64_t, ndim=2] vertices):
    
        b = np.empty(3, dtype=np.float64)
        A = np.empty((3, 3), dtype=np.float64)
    
        b[0] = physical_x[0]
        b[1] = physical_x[1]
        b[2] = 1.0
    
        A[0][0] = vertices[0, 0]
        A[0][1] = vertices[1, 0]
        A[0][2] = vertices[2, 0]
    
        A[1][0] = vertices[0, 1]
        A[1][1] = vertices[1, 1]
        A[1][2] = vertices[2, 1]
    
        A[2][0] = 1.0
        A[2][1] = 1.0
        A[2][2] = 1.0
            
        c = np.linalg.solve(A, b)
        return c

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def sample_at_unit_point(self, double[:] coord, 
                             double[:] vals):
        return vals[0]*coord[0] + vals[1]*coord[1] + vals[2]*coord[2]

cdef class P1Sampler3D(ElementSampler):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    def map_real_to_unit(self, double[:] physical_x, double[:,:] vertices):
    
        b = np.empty(4, dtype=np.float64)
        A = np.empty((4, 4), dtype=np.float64)
    
        b[0] = physical_x[0]
        b[1] = physical_x[1]
        b[2] = physical_x[2]
        b[3] = 1.0
    
        A[0][0] = vertices[0, 0]
        A[0][1] = vertices[1, 0]
        A[0][2] = vertices[2, 0]
        A[0][3] = vertices[3, 0]
        
        A[1][0] = vertices[0, 1]
        A[1][1] = vertices[1, 1]
        A[1][2] = vertices[2, 1]
        A[1][3] = vertices[3, 1]
        
        A[2][0] = vertices[0, 2]
        A[2][1] = vertices[1, 2]
        A[2][2] = vertices[2, 2]
        A[2][3] = vertices[3, 2]

        A[3][0] = 1.0
        A[3][1] = 1.0
        A[3][2] = 1.0
        A[3][3] = 1.0
        
        c = np.linalg.solve(A, b)
    
        return c

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def sample_at_unit_point(self,
                             double[:] coord, 
                             double[:] vals):
        cdef double value = 0.0
        cdef int i
        for i in range(4):
            value += vals[i]*coord[i]
        return value

ctypedef void (*func_type)(double[:], double[:], double[:, :], double[:])
ctypedef void (*jac_type)(double[:, :], double[:], double[:, :], double[:])

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef inline void Q1Function2D(double[:] fx,
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
cdef inline void Q1Jacobian2D(double[:, :] A,
                              double[:] x, 
                              double[:, :] v, 
                              double[:] phys_x) nogil:
    
    cdef int i
    for i in range(2):
        A[i][0] = -(1-x[1])*v[0][i] + (1-x[1])*v[1][i] - \
                   (1+x[1])*v[2][i] + (1+x[1])*v[3][i]
        A[i][1] = -(1-x[0])*v[0][i] - (1+x[0])*v[1][i] + \
                   (1-x[0])*v[2][i] + (1+x[0])*v[3][i]


cdef class Q1Sampler2D(NonlinearSolveSampler):

    def __init__(self):
        super(Q1Sampler2D, self).__init__()
        self.dim = 2
        self.func = Q1Function2D
        self.jac = Q1Jacobian2D

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
cdef inline void Q1Function3D(double[:] fx,
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
cdef inline void Q1Jacobian3D(double[:, :] A,
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

cdef class NonlinearSolveSampler(ElementSampler):

    cdef int dim
    cdef np.float64_t tolerance
    cdef func_type func 
    cdef jac_type jac

    def __init__(self):
        self.tolerance = 1.0e-9

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def map_real_to_unit(self, 
                         np.ndarray[np.float64_t, ndim=1] physical_x,
                         np.ndarray[np.float64_t, ndim=2] vertices):
        x = np.zeros(self.dim, dtype=np.float64)
        cdef int iterations = 0
        fx = np.empty(self.dim, dtype=np.float64)
        A = np.empty((self.dim, self.dim), dtype=np.float64)
        Ainv = np.empty((self.dim, self.dim), dtype=np.float64)
        self.func(fx, x, vertices, physical_x)
        cdef np.float64_t err = np.max(abs(fx))
        while (err > self.tolerance and iterations < 100):
            self.jac(A, x, vertices, physical_x)
            Ainv = np.linalg.inv(A)
            x = x - np.dot(Ainv, fx)
            self.func(fx, x, vertices, physical_x)
            err = np.max(abs(fx))
            iterations += 1
        return x


cdef class Q1Sampler3D(NonlinearSolveSampler):

    def __init__(self):
        super(Q1Sampler3D, self).__init__()
        self.dim = 3
        self.func = Q1Function3D
        self.jac = Q1Jacobian3D

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
