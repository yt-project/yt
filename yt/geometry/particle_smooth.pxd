"""
Particle Deposition onto Octs




"""


cimport numpy as np

import numpy as np

cimport cython
from libc.math cimport sqrt
from libc.stdlib cimport free, malloc, qsort
from oct_container cimport Oct, OctreeContainer

from yt.utilities.lib.distance_queue cimport (
    DistanceQueue,
    Neighbor_compare,
    NeighborList,
    r2dist,
)
from yt.utilities.lib.fp_utils cimport *

from .particle_deposit cimport get_kernel_func, gind, kernel_func


cdef extern from "platform_dep.h":
    void *alloca(int)

cdef class ParticleSmoothOperation:
    # We assume each will allocate and define their own temporary storage
    cdef kernel_func sph_kernel
    cdef public object nvals
    cdef np.float64_t DW[3]
    cdef int nfields
    cdef int maxn
    cdef bint periodicity[3]
    # Note that we are preallocating here, so this is *not* threadsafe.
    cdef void (*pos_setup)(np.float64_t ipos[3], np.float64_t opos[3])
    cdef void neighbor_process(self, int dim[3], np.float64_t left_edge[3],
                               np.float64_t dds[3], np.float64_t[:,:] ppos,
                               np.float64_t **fields,
                               np.int64_t[:] doffs, np.int64_t **nind,
                               np.int64_t[:] pinds, np.int64_t[:] pcounts,
                               np.int64_t offset, np.float64_t **index_fields,
                               OctreeContainer octree, np.int64_t domain_id,
                               int *nsize, np.float64_t[:,:] oct_left_edges,
                               np.float64_t[:,:] oct_dds, DistanceQueue dq)
    cdef int neighbor_search(self, np.float64_t pos[3], OctreeContainer octree,
                             np.int64_t **nind, int *nsize,
                             np.int64_t nneighbors, np.int64_t domain_id,
                             Oct **oct = ?, int extra_layer = ?)
    cdef void neighbor_process_particle(self, np.float64_t cpos[3],
                               np.float64_t[:,:] ppos,
                               np.float64_t **fields,
                               np.int64_t[:] doffs, np.int64_t **nind,
                               np.int64_t[:] pinds, np.int64_t[:] pcounts,
                               np.int64_t offset,
                               np.float64_t **index_fields,
                               OctreeContainer octree, np.int64_t domain_id,
                               int *nsize, DistanceQueue dq)
    cdef void neighbor_find(self,
                            np.int64_t nneighbors,
                            np.int64_t *nind,
                            np.int64_t[:] doffs,
                            np.int64_t[:] pcounts,
                            np.int64_t[:] pinds,
                            np.float64_t[:,:] ppos,
                            np.float64_t cpos[3],
                            np.float64_t[:,:] oct_left_edges,
                            np.float64_t[:,:] oct_dds, DistanceQueue dq)
    cdef void process(self, np.int64_t offset, int i, int j, int k,
                      int dim[3], np.float64_t cpos[3], np.float64_t **fields,
                      np.float64_t **index_fields, DistanceQueue dq)
