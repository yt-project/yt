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
from scipy.optimize import fsolve


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
cdef np.ndarray[np.float64_t, ndim=1] Q1Function2D(np.ndarray[np.float64_t, ndim=1] x,
                                                   np.ndarray[np.float64_t, ndim=2] v,
                                                   np.ndarray[np.float64_t, ndim=1] phys_x):
    f1 = v[0][0]*(1-x[0])*(1-x[1]) + \
         v[1][0]*(1+x[0])*(1-x[1]) + \
         v[2][0]*(1-x[0])*(1+x[1]) + \
         v[3][0]*(1+x[0])*(1+x[1]) - 4.0*phys_x[0]
    f2 = v[0][1]*(1-x[0])*(1-x[1]) + \
         v[1][1]*(1+x[0])*(1-x[1]) + \
         v[2][1]*(1-x[0])*(1+x[1]) + \
         v[3][1]*(1+x[0])*(1+x[1]) - 4.0*phys_x[1]
    return np.array([f1, f2])

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.ndarray[np.float64_t, ndim=2] Q1Jacobian2D(np.ndarray[np.float64_t, ndim=1] x,
                                                   np.ndarray[np.float64_t, ndim=2] v,
                                                   np.ndarray[np.float64_t, ndim=1] phys_x):
    f11 = -(1-x[1])*v[0][0] + \
          (1-x[1])*v[1][0] - \
          (1+x[1])*v[2][0] + \
          (1+x[1])*v[3][0]
    f12 = -(1-x[0])*v[0][0] - \
          (1+x[0])*v[1][0] + \
          (1-x[0])*v[2][0] + \
          (1+x[0])*v[3][0]
    f21 = -(1-x[1])*v[0][1] + \
          (1-x[1])*v[1][1] - \
          (1+x[1])*v[2][1] + \
          (1+x[1])*v[3][1]
    f22 = -(1-x[0])*v[0][1] - \
          (1+x[0])*v[1][1] + \
          (1-x[0])*v[2][1] + \
          (1+x[0])*v[3][1]
    return np.array([[f11, f12], [f21, f22]])


cdef class Q1Sampler2D(ElementSampler):

    def map_real_to_unit(self, np.ndarray[np.float64_t, ndim=1] physical_coord, 
                         np.ndarray[np.float64_t, ndim=2] vertices):
    
        # initial guess for the Newton solve
        x0 = np.array([0.0, 0.0], dtype=np.float64)
        x = fsolve(Q1Function2D, x0, args=(vertices, physical_coord),
                   fprime=Q1Jacobian2D)
        return x

    def sample_at_unit_point(self, np.ndarray[np.float64_t, ndim=1] coord, 
                             np.ndarray[np.float64_t, ndim=1] vals):
        cdef np.float64_t x = vals[0]*(1.0 - coord[0])*(1.0 - coord[1]) + \
                              vals[1]*(1.0 + coord[0])*(1.0 - coord[1]) + \
                              vals[2]*(1.0 - coord[0])*(1.0 + coord[1]) + \
                              vals[3]*(1.0 + coord[0])*(1.0 + coord[1])
        return 0.25*x

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline np.ndarray[np.float64_t, ndim=1] Q1Function3D(np.ndarray[np.float64_t, ndim=1] x,
                                                          np.ndarray[np.float64_t, ndim=2] v,
                                                          np.ndarray[np.float64_t, ndim=1] phys_x):
    cdef np.float64_t f0 = v[0][0]*(1-x[0])*(1-x[1])*(1-x[2]) + \
                           v[1][0]*(1+x[0])*(1-x[1])*(1-x[2]) + \
                           v[2][0]*(1-x[0])*(1+x[1])*(1-x[2]) + \
                           v[3][0]*(1+x[0])*(1+x[1])*(1-x[2]) + \
                           v[4][0]*(1-x[0])*(1-x[1])*(1+x[2]) + \
                           v[5][0]*(1+x[0])*(1-x[1])*(1+x[2]) + \
                           v[6][0]*(1-x[0])*(1+x[1])*(1+x[2]) + \
                           v[7][0]*(1+x[0])*(1+x[1])*(1+x[2]) - 8.0*phys_x[0]
    cdef np.float64_t f1 = v[0][1]*(1-x[0])*(1-x[1])*(1-x[2]) + \
                           v[1][1]*(1+x[0])*(1-x[1])*(1-x[2]) + \
                           v[2][1]*(1-x[0])*(1+x[1])*(1-x[2]) + \
                           v[3][1]*(1+x[0])*(1+x[1])*(1-x[2]) + \
                           v[4][1]*(1-x[0])*(1-x[1])*(1+x[2]) + \
                           v[5][1]*(1+x[0])*(1-x[1])*(1+x[2]) + \
                           v[6][1]*(1-x[0])*(1+x[1])*(1+x[2]) + \
                           v[7][1]*(1+x[0])*(1+x[1])*(1+x[2]) - 8.0*phys_x[1]
    cdef np.float64_t f2 = v[0][2]*(1-x[0])*(1-x[1])*(1-x[2]) + \
                           v[1][2]*(1+x[0])*(1-x[1])*(1-x[2]) + \
                           v[2][2]*(1-x[0])*(1+x[1])*(1-x[2]) + \
                           v[3][2]*(1+x[0])*(1+x[1])*(1-x[2]) + \
                           v[4][2]*(1-x[0])*(1-x[1])*(1+x[2]) + \
                           v[5][2]*(1+x[0])*(1-x[1])*(1+x[2]) + \
                           v[6][2]*(1-x[0])*(1+x[1])*(1+x[2]) + \
                           v[7][2]*(1+x[0])*(1+x[1])*(1+x[2]) - 8.0*phys_x[2]

    A = np.empty(3, dtype=np.float64)
    A[0] = f0
    A[1] = f1
    A[2] = f2
    return A

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.ndarray[np.float64_t, ndim=1] Q1Jacobian3D(np.ndarray[np.float64_t, ndim=1] x,
                                                   np.ndarray[np.float64_t, ndim=2] v,
                                                   np.ndarray[np.float64_t, ndim=1] phys_x):
    f00 = -(1-x[1])*(1-x[2])*v[0][0] + (1-x[1])*(1-x[2])*v[1][0] - \
          (1+x[1])*(1-x[2])*v[2][0] + (1+x[1])*(1-x[2])*v[3][0] - \
          (1-x[1])*(1+x[2])*v[4][0] + (1-x[1])*(1+x[2])*v[5][0] - \
          (1+x[1])*(1+x[2])*v[6][0] + (1+x[1])*(1+x[2])*v[7][0]
    f01 = -(1-x[0])*(1-x[2])*v[0][0] - (1+x[0])*(1-x[2])*v[1][0] + \
          (1-x[0])*(1-x[2])*v[2][0] + (1+x[0])*(1-x[2])*v[3][0] - \
          (1-x[0])*(1+x[2])*v[4][0] - (1+x[0])*(1+x[2])*v[5][0] + \
          (1-x[0])*(1+x[2])*v[6][0] + (1+x[0])*(1+x[2])*v[7][0]
    f02 = -(1-x[0])*(1-x[1])*v[0][0] - (1+x[0])*(1-x[1])*v[1][0] - \
          (1-x[0])*(1+x[1])*v[2][0] - (1+x[0])*(1+x[1])*v[3][0] + \
          (1-x[0])*(1-x[1])*v[4][0] + (1+x[0])*(1-x[1])*v[5][0] + \
          (1-x[0])*(1+x[1])*v[6][0] + (1+x[0])*(1+x[1])*v[7][0]

    f10 = -(1-x[1])*(1-x[2])*v[0][1] + (1-x[1])*(1-x[2])*v[1][1] - \
          (1+x[1])*(1-x[2])*v[2][1] + (1+x[1])*(1-x[2])*v[3][1] - \
          (1-x[1])*(1+x[2])*v[4][1] + (1-x[1])*(1+x[2])*v[5][1] - \
          (1+x[1])*(1+x[2])*v[6][1] + (1+x[1])*(1+x[2])*v[7][1]
    f11 = -(1-x[0])*(1-x[2])*v[0][1] - (1+x[0])*(1-x[2])*v[1][1] + \
          (1-x[0])*(1-x[2])*v[2][1] + (1+x[0])*(1-x[2])*v[3][1] - \
          (1-x[0])*(1+x[2])*v[4][1] - (1+x[0])*(1+x[2])*v[5][1] + \
          (1-x[0])*(1+x[2])*v[6][1] + (1+x[0])*(1+x[2])*v[7][1]
    f12 = -(1-x[0])*(1-x[1])*v[0][1] - (1+x[0])*(1-x[1])*v[1][1] - \
          (1-x[0])*(1+x[1])*v[2][1] - (1+x[0])*(1+x[1])*v[3][1] + \
          (1-x[0])*(1-x[1])*v[4][1] + (1+x[0])*(1-x[1])*v[5][1] + \
          (1-x[0])*(1+x[1])*v[6][1] + (1+x[0])*(1+x[1])*v[7][1]
    
    f20 = -(1-x[1])*(1-x[2])*v[0][2] + (1-x[1])*(1-x[2])*v[1][2] - \
          (1+x[1])*(1-x[2])*v[2][2] + (1+x[1])*(1-x[2])*v[3][2] - \
          (1-x[1])*(1+x[2])*v[4][2] + (1-x[1])*(1+x[2])*v[5][2] - \
          (1+x[1])*(1+x[2])*v[6][2] + (1+x[1])*(1+x[2])*v[7][2]
    f21 = -(1-x[0])*(1-x[2])*v[0][2] - (1+x[0])*(1-x[2])*v[1][2] + \
          (1-x[0])*(1-x[2])*v[2][2] + (1+x[0])*(1-x[2])*v[3][2] - \
          (1-x[0])*(1+x[2])*v[4][2] - (1+x[0])*(1+x[2])*v[5][2] + \
          (1-x[0])*(1+x[2])*v[6][2] + (1+x[0])*(1+x[2])*v[7][2]
    f22 = -(1-x[0])*(1-x[1])*v[0][2] - (1+x[0])*(1-x[1])*v[1][2] - \
          (1-x[0])*(1+x[1])*v[2][2] - (1+x[0])*(1+x[1])*v[3][2] + \
          (1-x[0])*(1-x[1])*v[4][2] + (1+x[0])*(1-x[1])*v[5][2] + \
          (1-x[0])*(1+x[1])*v[6][2] + (1+x[0])*(1+x[1])*v[7][2]

    return np.array([[f00, f01, f02], [f10, f11, f12], [f20, f21, f22]])


cdef class Q1Sampler3D(ElementSampler):

    def map_real_to_unit(self, np.ndarray physical_coord, np.ndarray vertices):
        x0 = np.array([0.0, 0.0, 0.0])  # initial guess
        x = fsolve(Q1Function3D, x0, args=(vertices, physical_coord),
                   fprime=Q1Jacobian3D)
        return x

    def sample_at_unit_point(self, np.ndarray coord, np.ndarray vals):
        x = vals[0]*(1.0 - coord[0])*(1.0 - coord[1])*(1.0 - coord[2]) + \
            vals[1]*(1.0 + coord[0])*(1.0 - coord[1])*(1.0 - coord[2]) + \
            vals[2]*(1.0 - coord[0])*(1.0 + coord[1])*(1.0 - coord[2]) + \
            vals[3]*(1.0 + coord[0])*(1.0 + coord[1])*(1.0 - coord[2]) + \
            vals[4]*(1.0 - coord[0])*(1.0 - coord[1])*(1.0 + coord[2]) + \
            vals[5]*(1.0 + coord[0])*(1.0 - coord[1])*(1.0 + coord[2]) + \
            vals[6]*(1.0 - coord[0])*(1.0 + coord[1])*(1.0 + coord[2]) + \
            vals[7]*(1.0 + coord[0])*(1.0 + coord[1])*(1.0 + coord[2])
        return 0.125*x
