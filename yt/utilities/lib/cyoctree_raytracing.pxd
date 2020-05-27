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

cdef extern from "octree_raytracing.cpp":
    cdef cppclass RayInfo[T]:
        vector[T] keys
        vector[double] t

    cdef cppclass Octree3D[T]:
        Octree3D(int depth, double* size)
        Octree3D(int depth, double* LE, double* RE)
        void insert_node_no_ret(const int* ipos, const int lvl, T key)
        RayInfo[T]** cast_rays(const double* origins, const double* directions, const int Nrays)

cdef class CythonOctreeRayTracing:
    cdef Octree3D[int]* oct
    cdef int depth
