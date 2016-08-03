# This file contains auto-generated functions for sampling 
# inside finite element solutions for various mesh types.
 
cimport cython 
 

 
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

 
