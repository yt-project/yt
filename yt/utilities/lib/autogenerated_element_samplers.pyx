# This file contains auto-generated functions for sampling 
# inside finite element solutions for various mesh types.
 
cimport cython 
 

 
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) 
cdef void Q1Function2D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil: 
    fx[0] =  -phys_x[0] + 0.25*vertices[0]*(-x[0] + 1)*(-x[1] + 1) + 0.25*vertices[2]*(x[0] + 1)*(-x[1] + 1) + 0.25*vertices[4]*(x[0] + 1)*(x[1] + 1) + 0.25*vertices[6]*(-x[0] + 1)*(x[1] + 1) 
    fx[1] =  -phys_x[1] + 0.25*vertices[1]*(-x[0] + 1)*(-x[1] + 1) + 0.25*vertices[3]*(x[0] + 1)*(-x[1] + 1) + 0.25*vertices[5]*(x[0] + 1)*(x[1] + 1) + 0.25*vertices[7]*(-x[0] + 1)*(x[1] + 1) 

 
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) 
cdef void Q1Jacobian2D(double* A,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil: 
    A[0] =  -0.25*vertices[0]*(-x[1] + 1) + 0.25*vertices[2]*(-x[1] + 1) + 0.25*vertices[4]*(x[1] + 1) - 0.25*vertices[6]*(x[1] + 1) 
    A[1] =  -0.25*vertices[0]*(-x[0] + 1) - 0.25*vertices[2]*(x[0] + 1) + 0.25*vertices[4]*(x[0] + 1) + 0.25*vertices[6]*(-x[0] + 1) 
    A[2] =  -0.25*vertices[1]*(-x[1] + 1) + 0.25*vertices[3]*(-x[1] + 1) + 0.25*vertices[5]*(x[1] + 1) - 0.25*vertices[7]*(x[1] + 1) 
    A[3] =  -0.25*vertices[1]*(-x[0] + 1) - 0.25*vertices[3]*(x[0] + 1) + 0.25*vertices[5]*(x[0] + 1) + 0.25*vertices[7]*(-x[0] + 1) 

 
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) 
cdef void W1Function3D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil: 
    fx[0] =  -phys_x[0] + 0.5*vertices[0]*(-x[2] + 1)*(-x[0] - x[1] + 1) + 0.5*vertices[3]*x[0]*(-x[2] + 1) + 0.5*vertices[6]*x[1]*(-x[2] + 1) + 0.5*vertices[9]*(x[2] + 1)*(-x[0] - x[1] + 1) + 0.5*vertices[12]*x[0]*(x[2] + 1) + 0.5*vertices[15]*x[1]*(x[2] + 1) 
    fx[1] =  -phys_x[1] + 0.5*vertices[1]*(-x[2] + 1)*(-x[0] - x[1] + 1) + 0.5*vertices[4]*x[0]*(-x[2] + 1) + 0.5*vertices[7]*x[1]*(-x[2] + 1) + 0.5*vertices[10]*(x[2] + 1)*(-x[0] - x[1] + 1) + 0.5*vertices[13]*x[0]*(x[2] + 1) + 0.5*vertices[16]*x[1]*(x[2] + 1) 
    fx[2] =  -phys_x[2] + 0.5*vertices[2]*(-x[2] + 1)*(-x[0] - x[1] + 1) + 0.5*vertices[5]*x[0]*(-x[2] + 1) + 0.5*vertices[8]*x[1]*(-x[2] + 1) + 0.5*vertices[11]*(x[2] + 1)*(-x[0] - x[1] + 1) + 0.5*vertices[14]*x[0]*(x[2] + 1) + 0.5*vertices[17]*x[1]*(x[2] + 1) 

 
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) 
cdef void W1Jacobian3D(double* rcol,
                       double* scol,
                       double* tcol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil: 
    rcol[0] =  -0.5*vertices[0]*(-x[2] + 1) + 0.5*vertices[3]*(-x[2] + 1) - 0.5*vertices[9]*(x[2] + 1) + 0.5*vertices[12]*(x[2] + 1) 
    scol[0] =  -0.5*vertices[0]*(-x[2] + 1) + 0.5*vertices[6]*(-x[2] + 1) - 0.5*vertices[9]*(x[2] + 1) + 0.5*vertices[15]*(x[2] + 1) 
    tcol[0] =  -0.5*vertices[0]*(-x[0] - x[1] + 1) - 0.5*vertices[3]*x[0] - 0.5*vertices[6]*x[1] + 0.5*vertices[9]*(-x[0] - x[1] + 1) + 0.5*vertices[12]*x[0] + 0.5*vertices[15]*x[1] 
    rcol[1] =  -0.5*vertices[1]*(-x[2] + 1) + 0.5*vertices[4]*(-x[2] + 1) - 0.5*vertices[10]*(x[2] + 1) + 0.5*vertices[13]*(x[2] + 1) 
    scol[1] =  -0.5*vertices[1]*(-x[2] + 1) + 0.5*vertices[7]*(-x[2] + 1) - 0.5*vertices[10]*(x[2] + 1) + 0.5*vertices[16]*(x[2] + 1) 
    tcol[1] =  -0.5*vertices[1]*(-x[0] - x[1] + 1) - 0.5*vertices[4]*x[0] - 0.5*vertices[7]*x[1] + 0.5*vertices[10]*(-x[0] - x[1] + 1) + 0.5*vertices[13]*x[0] + 0.5*vertices[16]*x[1] 
    rcol[2] =  -0.5*vertices[2]*(-x[2] + 1) + 0.5*vertices[5]*(-x[2] + 1) - 0.5*vertices[11]*(x[2] + 1) + 0.5*vertices[14]*(x[2] + 1) 
    scol[2] =  -0.5*vertices[2]*(-x[2] + 1) + 0.5*vertices[8]*(-x[2] + 1) - 0.5*vertices[11]*(x[2] + 1) + 0.5*vertices[17]*(x[2] + 1) 
    tcol[2] =  -0.5*vertices[2]*(-x[0] - x[1] + 1) - 0.5*vertices[5]*x[0] - 0.5*vertices[8]*x[1] + 0.5*vertices[11]*(-x[0] - x[1] + 1) + 0.5*vertices[14]*x[0] + 0.5*vertices[17]*x[1] 

 
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) 
cdef void Q1Function3D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil: 
    fx[0] =  -phys_x[0] + 0.125*vertices[0]*(-x[0] + 1)*(-x[1] + 1)*(-x[2] + 1) + 0.125*vertices[3]*(x[0] + 1)*(-x[1] + 1)*(-x[2] + 1) + 0.125*vertices[6]*(x[0] + 1)*(x[1] + 1)*(-x[2] + 1) + 0.125*vertices[9]*(-x[0] + 1)*(x[1] + 1)*(-x[2] + 1) + 0.125*vertices[12]*(-x[0] + 1)*(-x[1] + 1)*(x[2] + 1) + 0.125*vertices[15]*(x[0] + 1)*(-x[1] + 1)*(x[2] + 1) + 0.125*vertices[18]*(x[0] + 1)*(x[1] + 1)*(x[2] + 1) + 0.125*vertices[21]*(-x[0] + 1)*(x[1] + 1)*(x[2] + 1) 
    fx[1] =  -phys_x[1] + 0.125*vertices[1]*(-x[0] + 1)*(-x[1] + 1)*(-x[2] + 1) + 0.125*vertices[4]*(x[0] + 1)*(-x[1] + 1)*(-x[2] + 1) + 0.125*vertices[7]*(x[0] + 1)*(x[1] + 1)*(-x[2] + 1) + 0.125*vertices[10]*(-x[0] + 1)*(x[1] + 1)*(-x[2] + 1) + 0.125*vertices[13]*(-x[0] + 1)*(-x[1] + 1)*(x[2] + 1) + 0.125*vertices[16]*(x[0] + 1)*(-x[1] + 1)*(x[2] + 1) + 0.125*vertices[19]*(x[0] + 1)*(x[1] + 1)*(x[2] + 1) + 0.125*vertices[22]*(-x[0] + 1)*(x[1] + 1)*(x[2] + 1) 
    fx[2] =  -phys_x[2] + 0.125*vertices[2]*(-x[0] + 1)*(-x[1] + 1)*(-x[2] + 1) + 0.125*vertices[5]*(x[0] + 1)*(-x[1] + 1)*(-x[2] + 1) + 0.125*vertices[8]*(x[0] + 1)*(x[1] + 1)*(-x[2] + 1) + 0.125*vertices[11]*(-x[0] + 1)*(x[1] + 1)*(-x[2] + 1) + 0.125*vertices[14]*(-x[0] + 1)*(-x[1] + 1)*(x[2] + 1) + 0.125*vertices[17]*(x[0] + 1)*(-x[1] + 1)*(x[2] + 1) + 0.125*vertices[20]*(x[0] + 1)*(x[1] + 1)*(x[2] + 1) + 0.125*vertices[23]*(-x[0] + 1)*(x[1] + 1)*(x[2] + 1) 

 
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) 
cdef void Q1Jacobian3D(double* rcol,
                       double* scol,
                       double* tcol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil: 
    rcol[0] =  -0.125*vertices[0]*(-x[1] + 1)*(-x[2] + 1) + 0.125*vertices[3]*(-x[1] + 1)*(-x[2] + 1) + 0.125*vertices[6]*(x[1] + 1)*(-x[2] + 1) - 0.125*vertices[9]*(x[1] + 1)*(-x[2] + 1) - 0.125*vertices[12]*(-x[1] + 1)*(x[2] + 1) + 0.125*vertices[15]*(-x[1] + 1)*(x[2] + 1) + 0.125*vertices[18]*(x[1] + 1)*(x[2] + 1) - 0.125*vertices[21]*(x[1] + 1)*(x[2] + 1) 
    scol[0] =  -0.125*vertices[0]*(-x[0] + 1)*(-x[2] + 1) - 0.125*vertices[3]*(x[0] + 1)*(-x[2] + 1) + 0.125*vertices[6]*(x[0] + 1)*(-x[2] + 1) + 0.125*vertices[9]*(-x[0] + 1)*(-x[2] + 1) - 0.125*vertices[12]*(-x[0] + 1)*(x[2] + 1) - 0.125*vertices[15]*(x[0] + 1)*(x[2] + 1) + 0.125*vertices[18]*(x[0] + 1)*(x[2] + 1) + 0.125*vertices[21]*(-x[0] + 1)*(x[2] + 1) 
    tcol[0] =  -0.125*vertices[0]*(-x[0] + 1)*(-x[1] + 1) - 0.125*vertices[3]*(x[0] + 1)*(-x[1] + 1) - 0.125*vertices[6]*(x[0] + 1)*(x[1] + 1) - 0.125*vertices[9]*(-x[0] + 1)*(x[1] + 1) + 0.125*vertices[12]*(-x[0] + 1)*(-x[1] + 1) + 0.125*vertices[15]*(x[0] + 1)*(-x[1] + 1) + 0.125*vertices[18]*(x[0] + 1)*(x[1] + 1) + 0.125*vertices[21]*(-x[0] + 1)*(x[1] + 1) 
    rcol[1] =  -0.125*vertices[1]*(-x[1] + 1)*(-x[2] + 1) + 0.125*vertices[4]*(-x[1] + 1)*(-x[2] + 1) + 0.125*vertices[7]*(x[1] + 1)*(-x[2] + 1) - 0.125*vertices[10]*(x[1] + 1)*(-x[2] + 1) - 0.125*vertices[13]*(-x[1] + 1)*(x[2] + 1) + 0.125*vertices[16]*(-x[1] + 1)*(x[2] + 1) + 0.125*vertices[19]*(x[1] + 1)*(x[2] + 1) - 0.125*vertices[22]*(x[1] + 1)*(x[2] + 1) 
    scol[1] =  -0.125*vertices[1]*(-x[0] + 1)*(-x[2] + 1) - 0.125*vertices[4]*(x[0] + 1)*(-x[2] + 1) + 0.125*vertices[7]*(x[0] + 1)*(-x[2] + 1) + 0.125*vertices[10]*(-x[0] + 1)*(-x[2] + 1) - 0.125*vertices[13]*(-x[0] + 1)*(x[2] + 1) - 0.125*vertices[16]*(x[0] + 1)*(x[2] + 1) + 0.125*vertices[19]*(x[0] + 1)*(x[2] + 1) + 0.125*vertices[22]*(-x[0] + 1)*(x[2] + 1) 
    tcol[1] =  -0.125*vertices[1]*(-x[0] + 1)*(-x[1] + 1) - 0.125*vertices[4]*(x[0] + 1)*(-x[1] + 1) - 0.125*vertices[7]*(x[0] + 1)*(x[1] + 1) - 0.125*vertices[10]*(-x[0] + 1)*(x[1] + 1) + 0.125*vertices[13]*(-x[0] + 1)*(-x[1] + 1) + 0.125*vertices[16]*(x[0] + 1)*(-x[1] + 1) + 0.125*vertices[19]*(x[0] + 1)*(x[1] + 1) + 0.125*vertices[22]*(-x[0] + 1)*(x[1] + 1) 
    rcol[2] =  -0.125*vertices[2]*(-x[1] + 1)*(-x[2] + 1) + 0.125*vertices[5]*(-x[1] + 1)*(-x[2] + 1) + 0.125*vertices[8]*(x[1] + 1)*(-x[2] + 1) - 0.125*vertices[11]*(x[1] + 1)*(-x[2] + 1) - 0.125*vertices[14]*(-x[1] + 1)*(x[2] + 1) + 0.125*vertices[17]*(-x[1] + 1)*(x[2] + 1) + 0.125*vertices[20]*(x[1] + 1)*(x[2] + 1) - 0.125*vertices[23]*(x[1] + 1)*(x[2] + 1) 
    scol[2] =  -0.125*vertices[2]*(-x[0] + 1)*(-x[2] + 1) - 0.125*vertices[5]*(x[0] + 1)*(-x[2] + 1) + 0.125*vertices[8]*(x[0] + 1)*(-x[2] + 1) + 0.125*vertices[11]*(-x[0] + 1)*(-x[2] + 1) - 0.125*vertices[14]*(-x[0] + 1)*(x[2] + 1) - 0.125*vertices[17]*(x[0] + 1)*(x[2] + 1) + 0.125*vertices[20]*(x[0] + 1)*(x[2] + 1) + 0.125*vertices[23]*(-x[0] + 1)*(x[2] + 1) 
    tcol[2] =  -0.125*vertices[2]*(-x[0] + 1)*(-x[1] + 1) - 0.125*vertices[5]*(x[0] + 1)*(-x[1] + 1) - 0.125*vertices[8]*(x[0] + 1)*(x[1] + 1) - 0.125*vertices[11]*(-x[0] + 1)*(x[1] + 1) + 0.125*vertices[14]*(-x[0] + 1)*(-x[1] + 1) + 0.125*vertices[17]*(x[0] + 1)*(-x[1] + 1) + 0.125*vertices[20]*(x[0] + 1)*(x[1] + 1) + 0.125*vertices[23]*(-x[0] + 1)*(x[1] + 1) 

 
