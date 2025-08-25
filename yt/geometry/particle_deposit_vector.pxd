"""
Particle Deposition onto Octs




"""


cimport numpy as np

import numpy as np

cimport cython
from libc.math cimport sqrt
from libc.stdlib cimport free, malloc

from yt.utilities.lib.fp_utils cimport *

from .oct_container cimport Oct, OctreeContainer

cdef extern from "numpy/npy_math.h":
    double NPY_PI

cdef extern from "platform_dep.h":
    void *alloca(int)

cdef inline int gind(int i, int j, int k, int dims[3]):
    # The ordering is such that we want i to vary the slowest in this instance,
    # even though in other instances it varies the fastest.  To see this in
    # action, try looking at the results of an n_ref=256 particle CIC plot,
    # which shows it the most clearly.
    return ((i*dims[1])+j)*dims[2]+k


cdef class ParticleDepositOperation:
    # We assume each will allocate and define their own temporary storage
    cdef public object nvals
    cdef public int update_values
    cdef int process(self, int dim[3], int ipart, np.float64_t left_edge[3],
                     np.float64_t dds[3], np.int64_t offset,
                     np.float64_t ppos[3], np.float64_t[:, :] fields,
                     np.int64_t domain_ind) except -1 nogil
