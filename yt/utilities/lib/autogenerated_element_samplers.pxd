cdef void Q1Function3D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil


cdef void Q1Jacobian3D(double* rcol,
                       double* scol,
                       double* tcol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil


cdef void Q1Function2D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil


cdef void Q1Jacobian2D(double* rcol,
                       double* scol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil


cdef void Q2Function2D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil


cdef void Q2Jacobian2D(double* rcol,
                       double* scol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil


cdef void Tet2Function3D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil


cdef void Tet2Jacobian3D(double* rcol,
                       double* scol,
                       double* tcol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil


cdef void T2Function2D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil


cdef void T2Jacobian2D(double* rcol,
                       double* scol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil


cdef void W1Function3D(double* fx,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil


cdef void W1Jacobian3D(double* rcol,
                       double* scol,
                       double* tcol,
                       double* x,
                       double* vertices,
                       double* phys_x) nogil
