"""
Definitions for the partitioned grid




"""


import numpy as np

cimport cython
cimport numpy as np

from .volume_container cimport VolumeContainer


cdef class PartitionedGrid:
    cdef public object my_data
    cdef public object source_mask
    cdef public object LeftEdge
    cdef public object RightEdge
    cdef public int parent_grid_id
    cdef VolumeContainer *container
    cdef np.float64_t star_er
    cdef np.float64_t star_sigma_num
    cdef np.float64_t star_coeff
    cdef void get_vector_field(self, np.float64_t pos[3],
                               np.float64_t *vel, np.float64_t *vel_mag)
