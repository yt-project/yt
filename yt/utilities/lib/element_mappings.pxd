cimport numpy as np
from numpy cimport ndarray
cimport cython
import numpy as np
from libc.math cimport fabs, fmax

cdef class ElementSampler:

    # how close a point has to be to the element
    # to get counted as "inside". This is in the
    # mapped coordinates of the element.
    cdef np.float64_t inclusion_tol

    cdef void map_real_to_unit(self,
                               double* mapped_x, 
                               double* vertices,
                               double* physical_x) nogil


    cdef double sample_at_unit_point(self,
                                     double* coord,
                                     double* vals) nogil
    

    cdef double sample_at_real_point(self,
                                     double* vertices,
                                     double* field_values,
                                     double* physical_x) nogil

    cdef int check_inside(self, double* mapped_coord) nogil


cdef class P1Sampler3D(ElementSampler):

    cdef void map_real_to_unit(self,
                               double* mapped_x, 
                               double* vertices,
                               double* physical_x) nogil


    cdef double sample_at_unit_point(self,
                                     double* coord,
                                     double* vals) nogil

    cdef int check_inside(self, double* mapped_coord) nogil


ctypedef void (*func_type)(double*, double*, double*, double*) nogil
ctypedef void (*jac_type)(double*, double*, double*, double*, double*, double*) nogil

cdef class NonlinearSolveSampler(ElementSampler):

    cdef int dim
    cdef int max_iter
    cdef np.float64_t tolerance
    cdef func_type func 
    cdef jac_type jac

    cdef void map_real_to_unit(self,
                               double* mapped_x, 
                               double* vertices,
                               double* physical_x) nogil
    

cdef class Q1Sampler3D(NonlinearSolveSampler):

    cdef void map_real_to_unit(self,
                               double* mapped_x, 
                               double* vertices,
                               double* physical_x) nogil


    cdef double sample_at_unit_point(self,
                                     double* coord,
                                     double* vals) nogil

    cdef int check_inside(self, double* mapped_coord) nogil
