"""This is a wrapper around the C++ class to efficiently cast rays into an octree.
It relies on the seminal paper by  J. Revelles,, C.Ureña and M.Lastra.
"""


cimport numpy as np
import numpy as np
from libcpp.vector cimport vector
cimport cython
from cython.parallel import prange, parallel
from libc.stdlib cimport free, malloc

from .image_samplers cimport ImageAccumulator, ImageSampler
from .grid_traversal cimport sampler_function
from .volume_container cimport VolumeContainer
from .partitioned_grid cimport PartitionedGrid

DEF Nch = 4


cdef class CythonOctreeRayTracing:
    def __init__(self, np.ndarray LE, np.ndarray RE, int depth):
        cdef double* LE_ptr = <double *>LE.data
        cdef double* RE_ptr = <double *>RE.data
        self.oct = new Octree3D[int](depth, LE_ptr, RE_ptr)
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