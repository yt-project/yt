"""
Particle Deposition onto Octs




"""


cimport numpy as np

import numpy as np

cimport cython
from libc.math cimport sqrt
from libc.stdlib cimport free, malloc
from numpy.math cimport PI as NPY_PI

from yt.utilities.lib.fp_utils cimport *

from .oct_container cimport Oct, OctreeContainer


cdef extern from "platform_dep.h":
    void *alloca(int)

cdef inline int gind(int i, int j, int k, int dims[3]):
    # The ordering is such that we want i to vary the slowest in this instance,
    # even though in other instances it varies the fastest.  To see this in
    # action, try looking at the results of an n_ref=256 particle CIC plot,
    # which shows it the most clearly.
    return ((i*dims[1])+j)*dims[2]+k


####################################################
# Standard SPH kernel for use with the Grid method #
####################################################

cdef inline np.float64_t sph_kernel_cubic(np.float64_t x) nogil:
    cdef np.float64_t kernel
    # C is 8/pi
    cdef np.float64_t C = 2.5464790894703255
    if x <= 0.5:
        kernel = 1.-6.*x*x*(1.-x)
    elif x>0.5 and x<=1.0:
        kernel = 2.*(1.-x)*(1.-x)*(1.-x)
    else:
        kernel = 0.
    return kernel * C

########################################################
# Alternative SPH kernels for use with the Grid method #
########################################################

# quartic spline
cdef inline np.float64_t sph_kernel_quartic(np.float64_t x) nogil:
    cdef np.float64_t kernel
    cdef np.float64_t C = 9.71404681957369  # 5.**6/512/np.pi
    if x < 1:
        kernel = (1.-x)**4
        if x < 3./5:
            kernel -= 5*(3./5-x)**4
            if x < 1./5:
                kernel += 10*(1./5-x)**4
    else:
        kernel = 0.
    return kernel * C

# quintic spline
cdef inline np.float64_t sph_kernel_quintic(np.float64_t x) nogil:
    cdef np.float64_t kernel
    cdef np.float64_t C = 17.403593027098754  # 3.**7/40/np.pi
    if x < 1:
        kernel = (1.-x)**5
        if x < 2./3:
            kernel -= 6*(2./3-x)**5
            if x < 1./3:
                kernel += 15*(1./3-x)**5
    else:
        kernel = 0.
    return kernel * C

# Wendland C2
cdef inline np.float64_t sph_kernel_wendland2(np.float64_t x) nogil:
    cdef np.float64_t kernel
    cdef np.float64_t C = 3.3422538049298023  # 21./2/np.pi
    if x < 1:
        kernel = (1.-x)**4 * (1+4*x)
    else:
        kernel = 0.
    return kernel * C

# Wendland C4
cdef inline np.float64_t sph_kernel_wendland4(np.float64_t x) nogil:
    cdef np.float64_t kernel
    cdef np.float64_t C = 4.923856051905513  # 495./32/np.pi
    if x < 1:
        kernel = (1.-x)**6 * (1+6*x+35./3*x**2)
    else:
        kernel = 0.
    return kernel * C

# Wendland C6
cdef inline np.float64_t sph_kernel_wendland6(np.float64_t x) nogil:
    cdef np.float64_t kernel
    cdef np.float64_t C = 6.78895304126366  # 1365./64/np.pi
    if x < 1:
        kernel = (1.-x)**8 * (1+8*x+25*x**2+32*x**3)
    else:
        kernel = 0.
    return kernel * C

cdef inline np.float64_t sph_kernel_dummy(np.float64_t x) nogil:
    return 0

# I don't know the way to use a dict in a cdef class.
# So in order to mimic a registry functionality,
# I manually created a function to lookup the kernel functions.
ctypedef np.float64_t (*kernel_func) (np.float64_t) nogil
cdef inline kernel_func get_kernel_func(str kernel_name) nogil:
    with gil:
        if kernel_name == 'cubic':
            return sph_kernel_cubic
        elif kernel_name == 'quartic':
            return sph_kernel_quartic
        elif kernel_name == 'quintic':
            return sph_kernel_quintic
        elif kernel_name == 'wendland2':
            return sph_kernel_wendland2
        elif kernel_name == 'wendland4':
            return sph_kernel_wendland4
        elif kernel_name == 'wendland6':
            return sph_kernel_wendland6
        elif kernel_name == 'none':
            return sph_kernel_dummy
        else:
            raise NotImplementedError

cdef class ParticleDepositOperation:
    # We assume each will allocate and define their own temporary storage
    cdef kernel_func sph_kernel
    cdef public object nvals
    cdef public int update_values
    cdef int process(self, int dim[3], int ipart, np.float64_t left_edge[3],
                     np.float64_t dds[3], np.int64_t offset,
                     np.float64_t ppos[3], np.float64_t[:] fields,
                     np.int64_t domain_ind) nogil except -1
