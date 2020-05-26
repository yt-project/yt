"""This is a wrapper around the C++ class to efficiently cast rays into an octree.
It relies on the seminal paper by  J. Revelles,, C.Ureña and M.Lastra.
"""


cimport numpy as np
import numpy as np
from libcpp.vector cimport vector
cimport cython
from cython.parallel import prange, parallel
from libc.stdlib cimport free, malloc

from .image_samplers cimport ImageSampler, ImageAccumulator
from .volume_container cimport VolumeContainer

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

    @cython.boundscheck(False)
    @cython.wraparound(False)       
    def cast_rays(self, double[:, ::1] o, double[:, ::1] d, ImageSampler sampler, int num_threads = 0):
        cdef RayInfo[int]** ret
        cdef int Nrays = len(o)
        cdef RayInfo[int]* ri
        
        if Nrays == 0:
            return
        
        # print('Casting rays')
        
        ret = self.oct.cast_rays(&o[0,0], &d[0,0], Nrays)

        cdef int* cell_ind
        cdef double* tval

        cdef int i, j, vi, vj, nx, ny
        cdef VolumeContainer *vc
        cdef ImageAccumulator *idata
        cdef int* key_ptr
        cdef double* t_ptr
        cdef int[3] index = [1, 1, 1]

        nx = <int>np.round(Nrays**0.5)
        ny = nx  # TODO: change this

        with nogil, parallel(num_threads=num_threads):
            idata = <ImageAccumulator *> malloc(sizeof(ImageAccumulator))
            vc = <VolumeContainer*> malloc(sizeof(VolumeContainer))
            for j in prange(Nrays, schedule='static'):
                vj = j % ny
                vi = (j - vj) / ny
                ri = ret[j]
                if ri.keys.size() == 0:
                    continue

                for i in range(3):  # TODO: change 3 to Nchannel
                    idata.rgba[i] = 0

                key_ptr = &ri.keys[0]
                t_ptr = &ri.t[0]

                # Iterate over cells
                for i in range(ri.keys.size()):
                    # Now call the sampler on the list of cells
                    sampler.sample(
                        vc, 
                        &o[j, 0],
                        &d[j, 0],
                        t_ptr[2*i  ],
                        t_ptr[2*i+1],
                        index,
                        <void *> idata
                        )
        # Free memory
        for i in range(Nrays):
            free(ret[i])
        free(ret)

    def __dealloc__(self):
        del self.oct