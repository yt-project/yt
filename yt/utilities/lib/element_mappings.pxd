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
    cdef int num_mapped_coords

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

    cdef int check_mesh_lines(self, double* mapped_coord) nogil


cdef class P1Sampler2D(ElementSampler):

    cdef void map_real_to_unit(self,
                               double* mapped_x, 
                               double* vertices,
                               double* physical_x) nogil


    cdef double sample_at_unit_point(self,
                                     double* coord,
                                     double* vals) nogil

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

    cdef int check_mesh_lines(self, double* mapped_coord) nogil


# This typedef defines a function pointer that defines the system
# of equations that will be solved by the NonlinearSolveSamplers.
# 
# inputs:
#     x        - pointer to the mapped coordinate
#     vertices - pointer to the element vertices
#     phys_x   - pointer to the physical coordinate
#
# outputs:
#
#     fx - the result of solving the system, should be close to 0
#          once it is converged.
#
ctypedef void (*func_type)(double* fx, 
                           double* x, 
                           double* vertices, 
                           double* phys_x) nogil

# This typedef defines a function pointer that defines the Jacobian
# matrix used by the NonlinearSolveSampler3D. Subclasses needed to 
# define a Jacobian function in this form.
# 
# inputs:
#     x        - pointer to the mapped coordinate
#     vertices - pointer to the element vertices
#     phys_x   - pointer to the physical coordinate
#
# outputs:
#
#     rcol     - the first column of the jacobian
#     scol     - the second column of the jacobian
#     tcol     - the third column of the jaocobian
#
ctypedef void (*jac_type3D)(double* rcol, 
                            double* scol, 
                            double* tcol, 
                            double* x, 
                            double* vertices, 
                            double* phys_x) nogil


# This typedef defines a function pointer that defines the Jacobian
# matrix used by the NonlinearSolveSampler2D. Subclasses needed to 
# define a Jacobian function in this form.
# 
# inputs:
#     x        - pointer to the mapped coordinate
#     vertices - pointer to the element vertices
#     phys_x   - pointer to the physical coordinate
#
# outputs:
#
#     A        - A flattened array storing the Jacobian matrix
#                The order of this array is [J11, J12, J21, J22]
#
ctypedef void (*jac_type2D)(double* A,
                            double* x,
                            double* vertices,
                            double* phys_x) nogil


cdef class NonlinearSolveSampler3D(ElementSampler):

    cdef int dim
    cdef int max_iter
    cdef np.float64_t tolerance
    cdef func_type func 
    cdef jac_type3D jac

    cdef void map_real_to_unit(self,
                               double* mapped_x, 
                               double* vertices,
                               double* physical_x) nogil
    

cdef class Q1Sampler3D(NonlinearSolveSampler3D):

    cdef void map_real_to_unit(self,
                               double* mapped_x, 
                               double* vertices,
                               double* physical_x) nogil


    cdef double sample_at_unit_point(self,
                                     double* coord,
                                     double* vals) nogil

    cdef int check_inside(self, double* mapped_coord) nogil

    cdef int check_mesh_lines(self, double* mapped_coord) nogil


cdef class W1Sampler3D(NonlinearSolveSampler3D):

    cdef void map_real_to_unit(self,
                               double* mapped_x, 
                               double* vertices,
                               double* physical_x) nogil


    cdef double sample_at_unit_point(self,
                                     double* coord,
                                     double* vals) nogil

    cdef int check_inside(self, double* mapped_coord) nogil

    cdef int check_mesh_lines(self, double* mapped_coord) nogil



cdef class S2Sampler3D(NonlinearSolveSampler3D):

    cdef void map_real_to_unit(self,
                               double* mapped_x, 
                               double* vertices,
                               double* physical_x) nogil


    cdef double sample_at_unit_point(self,
                                     double* coord,
                                     double* vals) nogil

    cdef int check_inside(self, double* mapped_coord) nogil

    cdef int check_mesh_lines(self, double* mapped_coord) nogil



cdef class NonlinearSolveSampler2D(ElementSampler):

    cdef int dim
    cdef int max_iter
    cdef np.float64_t tolerance
    cdef func_type func 
    cdef jac_type2D jac

    cdef void map_real_to_unit(self,
                               double* mapped_x, 
                               double* vertices,
                               double* physical_x) nogil
    

cdef class Q1Sampler2D(NonlinearSolveSampler2D):

    cdef void map_real_to_unit(self,
                               double* mapped_x, 
                               double* vertices,
                               double* physical_x) nogil


    cdef double sample_at_unit_point(self,
                                     double* coord,
                                     double* vals) nogil

    cdef int check_inside(self, double* mapped_coord) nogil
