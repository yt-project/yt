"""
This file contains coordinate mappings between physical coordinates and those
defined on unit elements, as well as doing the corresponding intracell 
interpolation on finite element data.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np
from numpy cimport ndarray
cimport cython
import numpy as np
from libc.math cimport fabs, fmax


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double determinant_3x3(double* col0, 
                            double* col1, 
                            double* col2) nogil:
    return col0[0]*col1[1]*col2[2] - col0[0]*col1[2]*col2[1] - \
           col0[1]*col1[0]*col2[2] + col0[1]*col1[2]*col2[0] + \
           col0[2]*col1[0]*col2[1] - col0[2]*col1[1]*col2[0]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double maxnorm(double* f) nogil:
    cdef double err
    cdef int i
    err = fabs(f[0])
    for i in range(1, 2):
        err = fmax(err, fabs(f[i])) 
    return err


cdef class ElementSampler:
    '''

    This is a base class for sampling the value of a finite element solution
    at an arbitrary point inside a mesh element. In general, this will be done
    by transforming the requested physical coordinate into a mapped coordinate 
    system, sampling the solution in mapped coordinates, and returning the result.
    This is not to be used directly; use one of the subclasses instead.

    '''

    def __init__(self):
        self.inclusion_tol = 1.0e-8

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void map_real_to_unit(self,
                               double* mapped_x, 
                               double* vertices,
                               double* physical_x) nogil:
        pass
        
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef double sample_at_unit_point(self,
                                     double* coord,
                                     double* vals) nogil:
        pass

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int check_inside(self, double* mapped_coord) nogil:
        pass

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef double sample_at_real_point(self,
                                     double* vertices,
                                     double* field_values,
                                     double* physical_x) nogil:
        cdef double val
        cdef double mapped_coord[4]

        self.map_real_to_unit(mapped_coord, vertices, physical_x)
        val = self.sample_at_unit_point(mapped_coord, field_values)
        return val


cdef class P1Sampler3D(ElementSampler):
    '''

    This implements sampling inside a linear, tetrahedral mesh element.
    This mapping is linear and can be inverted easily.

    '''


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void map_real_to_unit(self, double* mapped_x, 
                               double* vertices, double* physical_x) nogil:
    
        cdef int i
        cdef double d
        cdef double[3] bvec
        cdef double[3] col0
        cdef double[3] col1
        cdef double[3] col2
    
        # here, we express positions relative to the 4th element,
        # which is selected by vertices[9]
        for i in range(3):
            bvec[i] = physical_x[i]       - vertices[9 + i]
            col0[i] = vertices[0 + i]     - vertices[9 + i]
            col1[i] = vertices[3 + i]     - vertices[9 + i]
            col2[i] = vertices[6 + i]     - vertices[9 + i]
        
        d = determinant_3x3(col0, col1, col2)
        mapped_x[0] = determinant_3x3(bvec, col1, col2)/d
        mapped_x[1] = determinant_3x3(col0, bvec, col2)/d
        mapped_x[2] = determinant_3x3(col0, col1, bvec)/d
        mapped_x[3] = 1.0 - mapped_x[0] - mapped_x[1] - mapped_x[2]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef double sample_at_unit_point(self,
                                     double* coord, 
                                     double* vals) nogil:
        return vals[0]*coord[0] + vals[1]*coord[1] + \
            vals[2]*coord[2] + vals[3]*coord[3]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int check_inside(self, double* mapped_coord) nogil:
        cdef int i
        for i in range(4):
            if (mapped_coord[i] < -self.inclusion_tol or
                mapped_coord[i] - 1.0 > self.inclusion_tol):
                return 0
        return 1


cdef class NonlinearSolveSampler(ElementSampler):

    '''

    This is a base class for handling element samplers that require
    a nonlinear solve to invert the mapping between coordinate systems.
    To do this, we perform Newton-Raphson iteration using a specified 
    system of equations with an analytic Jacobian matrix. This is
    not to be used directly, use one of the subclasses instead.

    '''

    def __init__(self):
        super(NonlinearSolveSampler, self).__init__()
        self.tolerance = 1.0e-9
        self.max_iter = 10

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void map_real_to_unit(self,
                               double* mapped_x,
                               double* vertices,
                               double* physical_x) nogil:
        cdef int i
        cdef double d, val
        cdef double[3] f
        cdef double[3] r
        cdef double[3] s
        cdef double[3] t
        cdef double[3] x
        cdef int iterations = 0
        cdef double err

        # initial guess
        for i in range(3):
            x[i] = 0.0
    
        # initial error norm
        self.func(f, x, vertices, physical_x)
        err = maxnorm(f)  
   
        # begin Newton iteration
        while (err > self.tolerance and iterations < self.max_iter):
            self.jac(r, s, t, x, vertices, physical_x)
            d = determinant_3x3(r, s, t)
            x[0] = x[0] - (determinant_3x3(f, s, t)/d)
            x[1] = x[1] - (determinant_3x3(r, f, t)/d)
            x[2] = x[2] - (determinant_3x3(r, s, f)/d)
            self.func(f, x, vertices, physical_x)        
            err = maxnorm(f)
            iterations += 1

        if (err > self.tolerance):
            # we did not converge, set bogus value
            for i in range(3):
                mapped_x[i] = -99.0
        else:
            for i in range(3):
                mapped_x[i] = x[i]


cdef class Q1Sampler3D(NonlinearSolveSampler):

    ''' 

    This implements sampling inside a 3D, linear, hexahedral mesh element.

    '''

    def __init__(self):
        super(Q1Sampler3D, self).__init__()
        self.dim = 3
        self.func = Q1Function3D
        self.jac = Q1Jacobian3D

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef double sample_at_unit_point(self, double* coord, double* vals) nogil:
        cdef double F, rm, rp, sm, sp, tm, tp
    
        rm = 1.0 - coord[0]
        rp = 1.0 + coord[0]
        sm = 1.0 - coord[1]
        sp = 1.0 + coord[1]
        tm = 1.0 - coord[2]
        tp = 1.0 + coord[2]
    
        F = vals[0]*rm*sm*tm + vals[1]*rp*sm*tm + vals[2]*rp*sp*tm + vals[3]*rm*sp*tm + \
            vals[4]*rm*sm*tp + vals[5]*rp*sm*tp + vals[6]*rp*sp*tp + vals[7]*rm*sp*tp
        return 0.125*F

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int check_inside(self, double* mapped_coord) nogil:
        if (fabs(mapped_coord[0]) - 1.0 > self.inclusion_tol or
            fabs(mapped_coord[1]) - 1.0 > self.inclusion_tol or 
            fabs(mapped_coord[2]) - 1.0 > self.inclusion_tol):
            return 0
        return 1


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void Q1Function3D(double* fx,
                              double* x, 
                              double* vertices, 
                              double* phys_x) nogil:
    cdef int i
    cdef double rm, rp, sm, sp, tm, tp
    
    rm = 1.0 - x[0]
    rp = 1.0 + x[0]
    sm = 1.0 - x[1]
    sp = 1.0 + x[1]
    tm = 1.0 - x[2]
    tp = 1.0 + x[2]
    
    for i in range(3):
        fx[i] = vertices[0 + i]*rm*sm*tm \
              + vertices[3 + i]*rp*sm*tm \
              + vertices[6 + i]*rp*sp*tm \
              + vertices[9 + i]*rm*sp*tm \
              + vertices[12 + i]*rm*sm*tp \
              + vertices[15 + i]*rp*sm*tp \
              + vertices[18 + i]*rp*sp*tp \
              + vertices[21 + i]*rm*sp*tp \
              - 8.0*phys_x[i]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void Q1Jacobian3D(double* rcol,
                              double* scol,
                              double* tcol,
                              double* x, 
                              double* vertices, 
                              double* phys_x) nogil:    
    cdef int i
    cdef double rm, rp, sm, sp, tm, tp
    
    rm = 1.0 - x[0]
    rp = 1.0 + x[0]
    sm = 1.0 - x[1]
    sp = 1.0 + x[1]
    tm = 1.0 - x[2]
    tp = 1.0 + x[2]
    
    for i in range(3):
        rcol[i] = -sm*tm*vertices[0 + i]  + sm*tm*vertices[3 + i]  + \
                   sp*tm*vertices[6 + i]  - sp*tm*vertices[9 + i]  - \
                   sm*tp*vertices[12 + i] + sm*tp*vertices[15 + i] + \
                   sp*tp*vertices[18 + i] - sp*tp*vertices[21 + i]
        scol[i] = -rm*tm*vertices[0 + i]  - rp*tm*vertices[3 + i]  + \
                   rp*tm*vertices[6 + i]  + rm*tm*vertices[9 + i]  - \
                   rm*tp*vertices[12 + i] - rp*tp*vertices[15 + i] + \
                   rp*tp*vertices[18 + i] + rm*tp*vertices[21 + i]
        tcol[i] = -rm*sm*vertices[0 + i]  - rp*sm*vertices[3 + i]  - \
                   rp*sp*vertices[6 + i]  - rm*sp*vertices[9 + i]  + \
                   rm*sm*vertices[12 + i] + rp*sm*vertices[15 + i] + \
                   rp*sp*vertices[18 + i] + rm*sp*vertices[21 + i]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def test_hex_sampler(np.ndarray[np.float64_t, ndim=2] vertices,
                     np.ndarray[np.float64_t, ndim=1] field_values,
                     np.ndarray[np.float64_t, ndim=1] physical_x):

    cdef double val

    cdef Q1Sampler3D sampler = Q1Sampler3D()

    val = sampler.sample_at_real_point(<double*> vertices.data,
                                       <double*> field_values.data,
                                       <double*> physical_x.data)
    return val


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def test_tetra_sampler(np.ndarray[np.float64_t, ndim=2] vertices,
                       np.ndarray[np.float64_t, ndim=1] field_values,
                       np.ndarray[np.float64_t, ndim=1] physical_x):

    cdef double val
    cdef double[4] mapped_coord

    sampler = P1Sampler3D()

    val = sampler.sample_at_real_point(<double*> vertices.data,
                                       <double*> field_values.data,
                                       <double*> physical_x.data)

    return val
