# This file contains auto-generated functions for sampling
# inside finite element solutions for various mesh types.
# To see how the code generation works in detail, see
# yt/utilities/mesh_code_generation.py.


cimport cython
from libc.math cimport pow


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void Q1Function3D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil:
	fx[0] = 0.125*(1 - x[0])*(1 - x[1])*(1 - x[2])*vertices[0] + 0.125*(1 - x[0])*(1 - x[1])*(1 + x[2])*vertices[12] + 0.125*(1 - x[0])*(1 + x[1])*(1 - x[2])*vertices[9] + 0.125*(1 - x[0])*(1 + x[1])*(1 + x[2])*vertices[21] + 0.125*(1 + x[0])*(1 - x[1])*(1 - x[2])*vertices[3] + 0.125*(1 + x[0])*(1 - x[1])*(1 + x[2])*vertices[15] + 0.125*(1 + x[0])*(1 + x[1])*(1 - x[2])*vertices[6] + 0.125*(1 + x[0])*(1 + x[1])*(1 + x[2])*vertices[18] - phys_x[0];
	fx[1] = 0.125*(1 - x[0])*(1 - x[1])*(1 - x[2])*vertices[1] + 0.125*(1 - x[0])*(1 - x[1])*(1 + x[2])*vertices[13] + 0.125*(1 - x[0])*(1 + x[1])*(1 - x[2])*vertices[10] + 0.125*(1 - x[0])*(1 + x[1])*(1 + x[2])*vertices[22] + 0.125*(1 + x[0])*(1 - x[1])*(1 - x[2])*vertices[4] + 0.125*(1 + x[0])*(1 - x[1])*(1 + x[2])*vertices[16] + 0.125*(1 + x[0])*(1 + x[1])*(1 - x[2])*vertices[7] + 0.125*(1 + x[0])*(1 + x[1])*(1 + x[2])*vertices[19] - phys_x[1];
	fx[2] = 0.125*(1 - x[0])*(1 - x[1])*(1 - x[2])*vertices[2] + 0.125*(1 - x[0])*(1 - x[1])*(1 + x[2])*vertices[14] + 0.125*(1 - x[0])*(1 + x[1])*(1 - x[2])*vertices[11] + 0.125*(1 - x[0])*(1 + x[1])*(1 + x[2])*vertices[23] + 0.125*(1 + x[0])*(1 - x[1])*(1 - x[2])*vertices[5] + 0.125*(1 + x[0])*(1 - x[1])*(1 + x[2])*vertices[17] + 0.125*(1 + x[0])*(1 + x[1])*(1 - x[2])*vertices[8] + 0.125*(1 + x[0])*(1 + x[1])*(1 + x[2])*vertices[20] - phys_x[2];


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void Q1Jacobian3D(double* rcol,
                       double* scol,
                       double* tcol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil:
	rcol[0] = -0.125*(1 - x[1])*(1 - x[2])*vertices[0] + 0.125*(1 - x[1])*(1 - x[2])*vertices[3] - 0.125*(1 - x[1])*(1 + x[2])*vertices[12] + 0.125*(1 - x[1])*(1 + x[2])*vertices[15] + 0.125*(1 + x[1])*(1 - x[2])*vertices[6] - 0.125*(1 + x[1])*(1 - x[2])*vertices[9] + 0.125*(1 + x[1])*(1 + x[2])*vertices[18] - 0.125*(1 + x[1])*(1 + x[2])*vertices[21];
	scol[0] = -0.125*(1 - x[0])*(1 - x[2])*vertices[0] + 0.125*(1 - x[0])*(1 - x[2])*vertices[9] - 0.125*(1 - x[0])*(1 + x[2])*vertices[12] + 0.125*(1 - x[0])*(1 + x[2])*vertices[21] - 0.125*(1 + x[0])*(1 - x[2])*vertices[3] + 0.125*(1 + x[0])*(1 - x[2])*vertices[6] - 0.125*(1 + x[0])*(1 + x[2])*vertices[15] + 0.125*(1 + x[0])*(1 + x[2])*vertices[18];
	tcol[0] = -0.125*(1 - x[0])*(1 - x[1])*vertices[0] + 0.125*(1 - x[0])*(1 - x[1])*vertices[12] - 0.125*(1 - x[0])*(1 + x[1])*vertices[9] + 0.125*(1 - x[0])*(1 + x[1])*vertices[21] - 0.125*(1 + x[0])*(1 - x[1])*vertices[3] + 0.125*(1 + x[0])*(1 - x[1])*vertices[15] - 0.125*(1 + x[0])*(1 + x[1])*vertices[6] + 0.125*(1 + x[0])*(1 + x[1])*vertices[18];
	rcol[1] = -0.125*(1 - x[1])*(1 - x[2])*vertices[1] + 0.125*(1 - x[1])*(1 - x[2])*vertices[4] - 0.125*(1 - x[1])*(1 + x[2])*vertices[13] + 0.125*(1 - x[1])*(1 + x[2])*vertices[16] + 0.125*(1 + x[1])*(1 - x[2])*vertices[7] - 0.125*(1 + x[1])*(1 - x[2])*vertices[10] + 0.125*(1 + x[1])*(1 + x[2])*vertices[19] - 0.125*(1 + x[1])*(1 + x[2])*vertices[22];
	scol[1] = -0.125*(1 - x[0])*(1 - x[2])*vertices[1] + 0.125*(1 - x[0])*(1 - x[2])*vertices[10] - 0.125*(1 - x[0])*(1 + x[2])*vertices[13] + 0.125*(1 - x[0])*(1 + x[2])*vertices[22] - 0.125*(1 + x[0])*(1 - x[2])*vertices[4] + 0.125*(1 + x[0])*(1 - x[2])*vertices[7] - 0.125*(1 + x[0])*(1 + x[2])*vertices[16] + 0.125*(1 + x[0])*(1 + x[2])*vertices[19];
	tcol[1] = -0.125*(1 - x[0])*(1 - x[1])*vertices[1] + 0.125*(1 - x[0])*(1 - x[1])*vertices[13] - 0.125*(1 - x[0])*(1 + x[1])*vertices[10] + 0.125*(1 - x[0])*(1 + x[1])*vertices[22] - 0.125*(1 + x[0])*(1 - x[1])*vertices[4] + 0.125*(1 + x[0])*(1 - x[1])*vertices[16] - 0.125*(1 + x[0])*(1 + x[1])*vertices[7] + 0.125*(1 + x[0])*(1 + x[1])*vertices[19];
	rcol[2] = -0.125*(1 - x[1])*(1 - x[2])*vertices[2] + 0.125*(1 - x[1])*(1 - x[2])*vertices[5] - 0.125*(1 - x[1])*(1 + x[2])*vertices[14] + 0.125*(1 - x[1])*(1 + x[2])*vertices[17] + 0.125*(1 + x[1])*(1 - x[2])*vertices[8] - 0.125*(1 + x[1])*(1 - x[2])*vertices[11] + 0.125*(1 + x[1])*(1 + x[2])*vertices[20] - 0.125*(1 + x[1])*(1 + x[2])*vertices[23];
	scol[2] = -0.125*(1 - x[0])*(1 - x[2])*vertices[2] + 0.125*(1 - x[0])*(1 - x[2])*vertices[11] - 0.125*(1 - x[0])*(1 + x[2])*vertices[14] + 0.125*(1 - x[0])*(1 + x[2])*vertices[23] - 0.125*(1 + x[0])*(1 - x[2])*vertices[5] + 0.125*(1 + x[0])*(1 - x[2])*vertices[8] - 0.125*(1 + x[0])*(1 + x[2])*vertices[17] + 0.125*(1 + x[0])*(1 + x[2])*vertices[20];
	tcol[2] = -0.125*(1 - x[0])*(1 - x[1])*vertices[2] + 0.125*(1 - x[0])*(1 - x[1])*vertices[14] - 0.125*(1 - x[0])*(1 + x[1])*vertices[11] + 0.125*(1 - x[0])*(1 + x[1])*vertices[23] - 0.125*(1 + x[0])*(1 - x[1])*vertices[5] + 0.125*(1 + x[0])*(1 - x[1])*vertices[17] - 0.125*(1 + x[0])*(1 + x[1])*vertices[8] + 0.125*(1 + x[0])*(1 + x[1])*vertices[20];


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void Q1Function2D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil:
	fx[0] = 0.25*(1 - x[0])*(1 - x[1])*vertices[0] + 0.25*(1 - x[0])*(1 + x[1])*vertices[6] + 0.25*(1 + x[0])*(1 - x[1])*vertices[2] + 0.25*(1 + x[0])*(1 + x[1])*vertices[4] - phys_x[0];
	fx[1] = 0.25*(1 - x[0])*(1 - x[1])*vertices[1] + 0.25*(1 - x[0])*(1 + x[1])*vertices[7] + 0.25*(1 + x[0])*(1 - x[1])*vertices[3] + 0.25*(1 + x[0])*(1 + x[1])*vertices[5] - phys_x[1];


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void Q1Jacobian2D(double* rcol,
                       double* scol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil:
	rcol[0] = -0.25*(1 - x[1])*vertices[0] + 0.25*(1 - x[1])*vertices[2] + 0.25*(1 + x[1])*vertices[4] - 0.25*(1 + x[1])*vertices[6];
	scol[0] = -0.25*(1 - x[0])*vertices[0] + 0.25*(1 - x[0])*vertices[6] - 0.25*(1 + x[0])*vertices[2] + 0.25*(1 + x[0])*vertices[4];
	rcol[1] = -0.25*(1 - x[1])*vertices[1] + 0.25*(1 - x[1])*vertices[3] + 0.25*(1 + x[1])*vertices[5] - 0.25*(1 + x[1])*vertices[7];
	scol[1] = -0.25*(1 - x[0])*vertices[1] + 0.25*(1 - x[0])*vertices[7] - 0.25*(1 + x[0])*vertices[3] + 0.25*(1 + x[0])*vertices[5];


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void Q2Function2D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil:
	fx[0] = (1 + x[0])*(-1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[16] - 0.5*(1 + x[0])*(-1 + x[0])*x[1]*(-1 + x[1])*vertices[8] - 0.5*(1 + x[0])*(-1 + x[0])*x[1]*(1 + x[1])*vertices[12] - phys_x[0] - 0.5*x[0]*(-1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[14] + 0.25*x[0]*(-1 + x[0])*x[1]*(-1 + x[1])*vertices[0] + 0.25*x[0]*(-1 + x[0])*x[1]*(1 + x[1])*vertices[6] - 0.5*x[0]*(1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[10] + 0.25*x[0]*(1 + x[0])*x[1]*(-1 + x[1])*vertices[2] + 0.25*x[0]*(1 + x[0])*x[1]*(1 + x[1])*vertices[4];
	fx[1] = (1 + x[0])*(-1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[17] - 0.5*(1 + x[0])*(-1 + x[0])*x[1]*(-1 + x[1])*vertices[9] - 0.5*(1 + x[0])*(-1 + x[0])*x[1]*(1 + x[1])*vertices[13] - phys_x[1] - 0.5*x[0]*(-1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[15] + 0.25*x[0]*(-1 + x[0])*x[1]*(-1 + x[1])*vertices[1] + 0.25*x[0]*(-1 + x[0])*x[1]*(1 + x[1])*vertices[7] - 0.5*x[0]*(1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[11] + 0.25*x[0]*(1 + x[0])*x[1]*(-1 + x[1])*vertices[3] + 0.25*x[0]*(1 + x[0])*x[1]*(1 + x[1])*vertices[5];


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void Q2Jacobian2D(double* rcol,
                       double* scol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil:
	rcol[0] = -0.5*(-1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[14] + (-1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[16] + 0.25*(-1 + x[0])*x[1]*(-1 + x[1])*vertices[0] - 0.5*(-1 + x[0])*x[1]*(-1 + x[1])*vertices[8] + 0.25*(-1 + x[0])*x[1]*(1 + x[1])*vertices[6] - 0.5*(-1 + x[0])*x[1]*(1 + x[1])*vertices[12] - 0.5*(1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[10] + (1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[16] + 0.25*(1 + x[0])*x[1]*(-1 + x[1])*vertices[2] - 0.5*(1 + x[0])*x[1]*(-1 + x[1])*vertices[8] + 0.25*(1 + x[0])*x[1]*(1 + x[1])*vertices[4] - 0.5*(1 + x[0])*x[1]*(1 + x[1])*vertices[12] - 0.5*x[0]*(1 + x[1])*(-1 + x[1])*vertices[10] - 0.5*x[0]*(1 + x[1])*(-1 + x[1])*vertices[14] + 0.25*x[0]*x[1]*(-1 + x[1])*vertices[0] + 0.25*x[0]*x[1]*(-1 + x[1])*vertices[2] + 0.25*x[0]*x[1]*(1 + x[1])*vertices[4] + 0.25*x[0]*x[1]*(1 + x[1])*vertices[6];
	scol[0] = -0.5*(-1 + x[1])*(1 + x[0])*(-1 + x[0])*vertices[8] + (-1 + x[1])*(1 + x[0])*(-1 + x[0])*vertices[16] + 0.25*(-1 + x[1])*x[0]*(-1 + x[0])*vertices[0] - 0.5*(-1 + x[1])*x[0]*(-1 + x[0])*vertices[14] + 0.25*(-1 + x[1])*x[0]*(1 + x[0])*vertices[2] - 0.5*(-1 + x[1])*x[0]*(1 + x[0])*vertices[10] - 0.5*(1 + x[1])*(1 + x[0])*(-1 + x[0])*vertices[12] + (1 + x[1])*(1 + x[0])*(-1 + x[0])*vertices[16] + 0.25*(1 + x[1])*x[0]*(-1 + x[0])*vertices[6] - 0.5*(1 + x[1])*x[0]*(-1 + x[0])*vertices[14] + 0.25*(1 + x[1])*x[0]*(1 + x[0])*vertices[4] - 0.5*(1 + x[1])*x[0]*(1 + x[0])*vertices[10] - 0.5*x[1]*(1 + x[0])*(-1 + x[0])*vertices[8] - 0.5*x[1]*(1 + x[0])*(-1 + x[0])*vertices[12] + 0.25*x[1]*x[0]*(-1 + x[0])*vertices[0] + 0.25*x[1]*x[0]*(-1 + x[0])*vertices[6] + 0.25*x[1]*x[0]*(1 + x[0])*vertices[2] + 0.25*x[1]*x[0]*(1 + x[0])*vertices[4];
	rcol[1] = -0.5*(-1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[15] + (-1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[17] + 0.25*(-1 + x[0])*x[1]*(-1 + x[1])*vertices[1] - 0.5*(-1 + x[0])*x[1]*(-1 + x[1])*vertices[9] + 0.25*(-1 + x[0])*x[1]*(1 + x[1])*vertices[7] - 0.5*(-1 + x[0])*x[1]*(1 + x[1])*vertices[13] - 0.5*(1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[11] + (1 + x[0])*(1 + x[1])*(-1 + x[1])*vertices[17] + 0.25*(1 + x[0])*x[1]*(-1 + x[1])*vertices[3] - 0.5*(1 + x[0])*x[1]*(-1 + x[1])*vertices[9] + 0.25*(1 + x[0])*x[1]*(1 + x[1])*vertices[5] - 0.5*(1 + x[0])*x[1]*(1 + x[1])*vertices[13] - 0.5*x[0]*(1 + x[1])*(-1 + x[1])*vertices[11] - 0.5*x[0]*(1 + x[1])*(-1 + x[1])*vertices[15] + 0.25*x[0]*x[1]*(-1 + x[1])*vertices[1] + 0.25*x[0]*x[1]*(-1 + x[1])*vertices[3] + 0.25*x[0]*x[1]*(1 + x[1])*vertices[5] + 0.25*x[0]*x[1]*(1 + x[1])*vertices[7];
	scol[1] = -0.5*(-1 + x[1])*(1 + x[0])*(-1 + x[0])*vertices[9] + (-1 + x[1])*(1 + x[0])*(-1 + x[0])*vertices[17] + 0.25*(-1 + x[1])*x[0]*(-1 + x[0])*vertices[1] - 0.5*(-1 + x[1])*x[0]*(-1 + x[0])*vertices[15] + 0.25*(-1 + x[1])*x[0]*(1 + x[0])*vertices[3] - 0.5*(-1 + x[1])*x[0]*(1 + x[0])*vertices[11] - 0.5*(1 + x[1])*(1 + x[0])*(-1 + x[0])*vertices[13] + (1 + x[1])*(1 + x[0])*(-1 + x[0])*vertices[17] + 0.25*(1 + x[1])*x[0]*(-1 + x[0])*vertices[7] - 0.5*(1 + x[1])*x[0]*(-1 + x[0])*vertices[15] + 0.25*(1 + x[1])*x[0]*(1 + x[0])*vertices[5] - 0.5*(1 + x[1])*x[0]*(1 + x[0])*vertices[11] - 0.5*x[1]*(1 + x[0])*(-1 + x[0])*vertices[9] - 0.5*x[1]*(1 + x[0])*(-1 + x[0])*vertices[13] + 0.25*x[1]*x[0]*(-1 + x[0])*vertices[1] + 0.25*x[1]*x[0]*(-1 + x[0])*vertices[7] + 0.25*x[1]*x[0]*(1 + x[0])*vertices[3] + 0.25*x[1]*x[0]*(1 + x[0])*vertices[5];


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void Tet2Function3D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil:
	fx[0] = (-x[0] + 2*pow(x[0], 2))*vertices[3] + (-x[1] + 2*pow(x[1], 2))*vertices[6] + (-x[2] + 2*pow(x[2], 2))*vertices[9] + (4*x[0] - 4*x[0]*x[1] - 4*x[0]*x[2] - 4*pow(x[0], 2))*vertices[12] + (4*x[1] - 4*x[1]*x[0] - 4*x[1]*x[2] - 4*pow(x[1], 2))*vertices[18] + (4*x[2] - 4*x[2]*x[0] - 4*x[2]*x[1] - 4*pow(x[2], 2))*vertices[21] + (1 - 3*x[0] + 4*x[0]*x[1] + 4*x[0]*x[2] + 2*pow(x[0], 2) - 3*x[1] + 4*x[1]*x[2] + 2*pow(x[1], 2) - 3*x[2] + 2*pow(x[2], 2))*vertices[0] - phys_x[0] + 4*x[0]*x[1]*vertices[15] + 4*x[0]*x[2]*vertices[24] + 4*x[1]*x[2]*vertices[27];
	fx[1] = (-x[0] + 2*pow(x[0], 2))*vertices[4] + (-x[1] + 2*pow(x[1], 2))*vertices[7] + (-x[2] + 2*pow(x[2], 2))*vertices[10] + (4*x[0] - 4*x[0]*x[1] - 4*x[0]*x[2] - 4*pow(x[0], 2))*vertices[13] + (4*x[1] - 4*x[1]*x[0] - 4*x[1]*x[2] - 4*pow(x[1], 2))*vertices[19] + (4*x[2] - 4*x[2]*x[0] - 4*x[2]*x[1] - 4*pow(x[2], 2))*vertices[22] + (1 - 3*x[0] + 4*x[0]*x[1] + 4*x[0]*x[2] + 2*pow(x[0], 2) - 3*x[1] + 4*x[1]*x[2] + 2*pow(x[1], 2) - 3*x[2] + 2*pow(x[2], 2))*vertices[1] - phys_x[1] + 4*x[0]*x[1]*vertices[16] + 4*x[0]*x[2]*vertices[25] + 4*x[1]*x[2]*vertices[28];
	fx[2] = (-x[0] + 2*pow(x[0], 2))*vertices[5] + (-x[1] + 2*pow(x[1], 2))*vertices[8] + (-x[2] + 2*pow(x[2], 2))*vertices[11] + (4*x[0] - 4*x[0]*x[1] - 4*x[0]*x[2] - 4*pow(x[0], 2))*vertices[14] + (4*x[1] - 4*x[1]*x[0] - 4*x[1]*x[2] - 4*pow(x[1], 2))*vertices[20] + (4*x[2] - 4*x[2]*x[0] - 4*x[2]*x[1] - 4*pow(x[2], 2))*vertices[23] + (1 - 3*x[0] + 4*x[0]*x[1] + 4*x[0]*x[2] + 2*pow(x[0], 2) - 3*x[1] + 4*x[1]*x[2] + 2*pow(x[1], 2) - 3*x[2] + 2*pow(x[2], 2))*vertices[2] - phys_x[2] + 4*x[0]*x[1]*vertices[17] + 4*x[0]*x[2]*vertices[26] + 4*x[1]*x[2]*vertices[29];


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void Tet2Jacobian3D(double* rcol,
                       double* scol,
                       double* tcol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil:
	rcol[0] = (-1 + 4*x[0])*vertices[3] + (-3 + 4*x[0] + 4*x[1] + 4*x[2])*vertices[0] + (4 - 8*x[0] - 4*x[1] - 4*x[2])*vertices[12] + 4*x[1]*vertices[15] - 4*x[1]*vertices[18] - 4*x[2]*vertices[21] + 4*x[2]*vertices[24];
	scol[0] = (-1 + 4*x[1])*vertices[6] + (-3 + 4*x[0] + 4*x[1] + 4*x[2])*vertices[0] + (4 - 4*x[0] - 8*x[1] - 4*x[2])*vertices[18] - 4*x[0]*vertices[12] + 4*x[0]*vertices[15] - 4*x[2]*vertices[21] + 4*x[2]*vertices[27];
	tcol[0] = (-1 + 4*x[2])*vertices[9] + (-3 + 4*x[0] + 4*x[1] + 4*x[2])*vertices[0] + (4 - 4*x[0] - 4*x[1] - 8*x[2])*vertices[21] - 4*x[0]*vertices[12] + 4*x[0]*vertices[24] - 4*x[1]*vertices[18] + 4*x[1]*vertices[27];
	rcol[1] = (-1 + 4*x[0])*vertices[4] + (-3 + 4*x[0] + 4*x[1] + 4*x[2])*vertices[1] + (4 - 8*x[0] - 4*x[1] - 4*x[2])*vertices[13] + 4*x[1]*vertices[16] - 4*x[1]*vertices[19] - 4*x[2]*vertices[22] + 4*x[2]*vertices[25];
	scol[1] = (-1 + 4*x[1])*vertices[7] + (-3 + 4*x[0] + 4*x[1] + 4*x[2])*vertices[1] + (4 - 4*x[0] - 8*x[1] - 4*x[2])*vertices[19] - 4*x[0]*vertices[13] + 4*x[0]*vertices[16] - 4*x[2]*vertices[22] + 4*x[2]*vertices[28];
	tcol[1] = (-1 + 4*x[2])*vertices[10] + (-3 + 4*x[0] + 4*x[1] + 4*x[2])*vertices[1] + (4 - 4*x[0] - 4*x[1] - 8*x[2])*vertices[22] - 4*x[0]*vertices[13] + 4*x[0]*vertices[25] - 4*x[1]*vertices[19] + 4*x[1]*vertices[28];
	rcol[2] = (-1 + 4*x[0])*vertices[5] + (-3 + 4*x[0] + 4*x[1] + 4*x[2])*vertices[2] + (4 - 8*x[0] - 4*x[1] - 4*x[2])*vertices[14] + 4*x[1]*vertices[17] - 4*x[1]*vertices[20] - 4*x[2]*vertices[23] + 4*x[2]*vertices[26];
	scol[2] = (-1 + 4*x[1])*vertices[8] + (-3 + 4*x[0] + 4*x[1] + 4*x[2])*vertices[2] + (4 - 4*x[0] - 8*x[1] - 4*x[2])*vertices[20] - 4*x[0]*vertices[14] + 4*x[0]*vertices[17] - 4*x[2]*vertices[23] + 4*x[2]*vertices[29];
	tcol[2] = (-1 + 4*x[2])*vertices[11] + (-3 + 4*x[0] + 4*x[1] + 4*x[2])*vertices[2] + (4 - 4*x[0] - 4*x[1] - 8*x[2])*vertices[23] - 4*x[0]*vertices[14] + 4*x[0]*vertices[26] - 4*x[1]*vertices[20] + 4*x[1]*vertices[29];


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void T2Function2D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil:
	fx[0] = (-x[0] + 2*pow(x[0], 2))*vertices[2] + (-x[1] + 2*pow(x[1], 2))*vertices[4] + (-4*x[0]*x[1] + 4*x[1] - 4*pow(x[1], 2))*vertices[10] + (4*x[0] - 4*x[0]*x[1] - 4*pow(x[0], 2))*vertices[6] + (1 - 3*x[0] + 4*x[0]*x[1] + 2*pow(x[0], 2) - 3*x[1] + 2*pow(x[1], 2))*vertices[0] - phys_x[0] + 4*x[0]*x[1]*vertices[8];
	fx[1] = (-x[0] + 2*pow(x[0], 2))*vertices[3] + (-x[1] + 2*pow(x[1], 2))*vertices[5] + (-4*x[0]*x[1] + 4*x[1] - 4*pow(x[1], 2))*vertices[11] + (4*x[0] - 4*x[0]*x[1] - 4*pow(x[0], 2))*vertices[7] + (1 - 3*x[0] + 4*x[0]*x[1] + 2*pow(x[0], 2) - 3*x[1] + 2*pow(x[1], 2))*vertices[1] - phys_x[1] + 4*x[0]*x[1]*vertices[9];


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void T2Jacobian2D(double* rcol,
                       double* scol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil:
	rcol[0] = (-1 + 4*x[0])*vertices[2] + (-3 + 4*x[0] + 4*x[1])*vertices[0] + (4 - 8*x[0] - 4*x[1])*vertices[6] + 4*x[1]*vertices[8] - 4*x[1]*vertices[10];
	scol[0] = (-1 + 4*x[1])*vertices[4] + (-3 + 4*x[0] + 4*x[1])*vertices[0] + (4 - 4*x[0] - 8*x[1])*vertices[10] - 4*x[0]*vertices[6] + 4*x[0]*vertices[8];
	rcol[1] = (-1 + 4*x[0])*vertices[3] + (-3 + 4*x[0] + 4*x[1])*vertices[1] + (4 - 8*x[0] - 4*x[1])*vertices[7] + 4*x[1]*vertices[9] - 4*x[1]*vertices[11];
	scol[1] = (-1 + 4*x[1])*vertices[5] + (-3 + 4*x[0] + 4*x[1])*vertices[1] + (4 - 4*x[0] - 8*x[1])*vertices[11] - 4*x[0]*vertices[7] + 4*x[0]*vertices[9];


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void W1Function3D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil:
	fx[0] = 0.5*(1 - x[0] - x[1])*(1 - x[2])*vertices[0] + 0.5*(1 - x[0] - x[1])*(1 + x[2])*vertices[9] - phys_x[0] + 0.5*x[0]*(1 - x[2])*vertices[3] + 0.5*x[0]*(1 + x[2])*vertices[12] + 0.5*x[1]*(1 - x[2])*vertices[6] + 0.5*x[1]*(1 + x[2])*vertices[15];
	fx[1] = 0.5*(1 - x[0] - x[1])*(1 - x[2])*vertices[1] + 0.5*(1 - x[0] - x[1])*(1 + x[2])*vertices[10] - phys_x[1] + 0.5*x[0]*(1 - x[2])*vertices[4] + 0.5*x[0]*(1 + x[2])*vertices[13] + 0.5*x[1]*(1 - x[2])*vertices[7] + 0.5*x[1]*(1 + x[2])*vertices[16];
	fx[2] = 0.5*(1 - x[0] - x[1])*(1 - x[2])*vertices[2] + 0.5*(1 - x[0] - x[1])*(1 + x[2])*vertices[11] - phys_x[2] + 0.5*x[0]*(1 - x[2])*vertices[5] + 0.5*x[0]*(1 + x[2])*vertices[14] + 0.5*x[1]*(1 - x[2])*vertices[8] + 0.5*x[1]*(1 + x[2])*vertices[17];


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void W1Jacobian3D(double* rcol,
                       double* scol,
                       double* tcol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil:
	rcol[0] = -0.5*(1 - x[2])*vertices[0] + 0.5*(1 - x[2])*vertices[3] - 0.5*(1 + x[2])*vertices[9] + 0.5*(1 + x[2])*vertices[12];
	scol[0] = -0.5*(1 - x[2])*vertices[0] + 0.5*(1 - x[2])*vertices[6] - 0.5*(1 + x[2])*vertices[9] + 0.5*(1 + x[2])*vertices[15];
	tcol[0] = -0.5*(1 - x[0] - x[1])*vertices[0] + 0.5*(1 - x[0] - x[1])*vertices[9] - 0.5*x[0]*vertices[3] + 0.5*x[0]*vertices[12] - 0.5*x[1]*vertices[6] + 0.5*x[1]*vertices[15];
	rcol[1] = -0.5*(1 - x[2])*vertices[1] + 0.5*(1 - x[2])*vertices[4] - 0.5*(1 + x[2])*vertices[10] + 0.5*(1 + x[2])*vertices[13];
	scol[1] = -0.5*(1 - x[2])*vertices[1] + 0.5*(1 - x[2])*vertices[7] - 0.5*(1 + x[2])*vertices[10] + 0.5*(1 + x[2])*vertices[16];
	tcol[1] = -0.5*(1 - x[0] - x[1])*vertices[1] + 0.5*(1 - x[0] - x[1])*vertices[10] - 0.5*x[0]*vertices[4] + 0.5*x[0]*vertices[13] - 0.5*x[1]*vertices[7] + 0.5*x[1]*vertices[16];
	rcol[2] = -0.5*(1 - x[2])*vertices[2] + 0.5*(1 - x[2])*vertices[5] - 0.5*(1 + x[2])*vertices[11] + 0.5*(1 + x[2])*vertices[14];
	scol[2] = -0.5*(1 - x[2])*vertices[2] + 0.5*(1 - x[2])*vertices[8] - 0.5*(1 + x[2])*vertices[11] + 0.5*(1 + x[2])*vertices[17];
	tcol[2] = -0.5*(1 - x[0] - x[1])*vertices[2] + 0.5*(1 - x[0] - x[1])*vertices[11] - 0.5*x[0]*vertices[5] + 0.5*x[0]*vertices[14] - 0.5*x[1]*vertices[8] + 0.5*x[1]*vertices[17];
