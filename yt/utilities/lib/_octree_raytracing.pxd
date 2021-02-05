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


cdef extern from "_octree_raytracing.hpp":
    cdef cppclass RayInfo[T]:
        vector[T] keys
        vector[double] t

    cdef cppclass Octree[T] nogil:
        Octree(int depth, double* LE, double* RE)
        void insert_node_no_ret(const int* ipos, const int lvl, T key)
        void cast_ray(double* origins, double* directions, vector[T] keyList, vector[double] tList)

cdef class _OctreeRayTracing:
    cdef Octree[int]* oct
    cdef int depth
