# distutils: language = c++
# distutils: extra_compile_args = CPP14_FLAG
"""This is a wrapper around the C++ class to efficiently cast rays into an octree.
It relies on the seminal paper by  J. Revelles,, C.Ureña and M.Lastra.
"""


cimport numpy as np

import numpy as np

cimport cython
from libcpp.vector cimport vector

from cython.parallel import parallel, prange

from libc.stdlib cimport free, malloc

from .grid_traversal cimport sampler_function
from .image_samplers cimport ImageAccumulator, ImageSampler
from .partitioned_grid cimport PartitionedGrid
from .volume_container cimport VolumeContainer

DEF Nch = 4


cdef class _OctreeRayTracing:
    def __init__(self, np.ndarray LE, np.ndarray RE, int depth):
        cdef double* LE_ptr = <double *>LE.data
        cdef double* RE_ptr = <double *>RE.data
        self.oct = new Octree[int](depth, LE_ptr, RE_ptr)
        self.depth = depth

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def add_nodes(self, int[:, :] ipos_view, int[:] lvl_view, int[:] key):
        cdef int i
        cdef int ii[3]

        for i in range(len(key)):
            ii[0] = ipos_view[i, 0]
            ii[1] = ipos_view[i, 1]
            ii[2] = ipos_view[i, 2]
            self.oct.insert_node_no_ret(ii, lvl_view[i], <int> key[i])

    def __dealloc__(self):
        del self.oct
