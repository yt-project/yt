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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def cast_rays(self, double[:, ::1] o, double[:, ::1] d, ImageSampler sampler, PartitionedGrid pg, int num_threads = 0):
        cdef RayInfo[int]** ret
        cdef int Nrays = len(o)
        cdef RayInfo[int]* ri

        cdef sampler_function* sample = <sampler_function *>sampler.sample

        if Nrays == 0:
            return

        ret = self.oct.cast_rays(&o[0,0], &d[0,0], Nrays)

        cdef int* cell_ind
        cdef double* tval

        cdef int i, j, k, vi, vj, nx, ny, icell
        cdef VolumeContainer *vc
        cdef ImageAccumulator *idata
        cdef int* key_ptr
        cdef double* t_ptr
        cdef int[3] index = [0, 0, 0]

        # Only support square regions for the moment!
        nx = <int>np.round(Nrays**0.5)
        ny = nx

        cdef int n_fields = pg.container.n_fields

        with nogil, parallel(num_threads=num_threads):
            idata = <ImageAccumulator *> malloc(sizeof(ImageAccumulator))
            vc = <VolumeContainer*> malloc(sizeof(VolumeContainer))
            vc.n_fields = 1
            vc.data = <np.float64_t**> malloc(sizeof(np.float64_t*))
            vc.mask = <np.uint8_t*> malloc(8*sizeof(np.uint8_t))
            # The actual dimensions are 2x2x2, but the sampler
            # assumes vertex-centred data for a 1x1x1 lattice (i.e.
            # 2^3 vertices)
            vc.dims[0] = 1
            vc.dims[1] = 1
            vc.dims[2] = 1
            for j in prange(Nrays, schedule='static'):
                vj = j % ny
                vi = (j - vj) / ny
                ri = ret[j]
                if ri.keys.size() == 0:
                    continue

                for i in range(Nch):
                    idata.rgba[i] = 0
                for i in range(8):
                    vc.mask[i] = 1
                key_ptr = &ri.keys[0]
                t_ptr = &ri.t[0]

                # Iterate over cells
                for i in range(ri.keys.size()):
                    icell = key_ptr[i]
                    for k in range(n_fields):
                        vc.data[k] = &pg.container.data[k][8*icell]

                    # Fill the volume container
                    for k in range(3):
                        vc.left_edge[k] = pg.container.left_edge[3*icell+k]
                        vc.right_edge[k] = pg.container.right_edge[3*icell+k]
                        vc.dds[i] = (vc.right_edge[i] - vc.left_edge[i])
                        vc.idds[i] = 1/vc.dds[i]
                    # Now call the sampler on the list of cells
                    sample(
                        vc,
                        &o[j, 0],
                        &d[j, 0],
                        t_ptr[2*i  ],
                        t_ptr[2*i+1],
                        index,
                        <void *> idata
                        )

                for i in range(Nch):
                    idata.rgba[i] = 0
            free(vc.data)
            free(vc)
            free(idata)
        # Free memory
        for i in range(Nrays):
            free(ret[i])
        free(ret)

    def __dealloc__(self):
        del self.oct