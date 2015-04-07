"""
Definitions for the traversal code




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython
cimport kdtree_utils

cdef struct VolumeContainer:
    int n_fields
    np.float64_t **data
    # The mask has dimensions one fewer in each direction than data
    np.uint8_t *mask
    np.float64_t left_edge[3]
    np.float64_t right_edge[3]
    np.float64_t dds[3]
    np.float64_t idds[3]
    int dims[3]

cdef class PartitionedGrid:
    cdef public object my_data
    cdef public object source_mask
    cdef public object LeftEdge
    cdef public object RightEdge
    cdef public int parent_grid_id
    cdef VolumeContainer *container
    cdef kdtree_utils.kdtree *star_list
    cdef np.float64_t star_er
    cdef np.float64_t star_sigma_num
    cdef np.float64_t star_coeff
    cdef void get_vector_field(self, np.float64_t pos[3],
                               np.float64_t *vel, np.float64_t *vel_mag)

ctypedef void sample_function(
                VolumeContainer *vc,
                np.float64_t v_pos[3],
                np.float64_t v_dir[3],
                np.float64_t enter_t,
                np.float64_t exit_t,
                int index[3],
                void *data) nogil

cdef int walk_volume(VolumeContainer *vc,
                     np.float64_t v_pos[3],
                     np.float64_t v_dir[3],
                     sample_function *sampler,
                     void *data,
                     np.float64_t *return_t = *,
                     np.float64_t enter_t = *) nogil

cdef inline int vc_index(VolumeContainer *vc, int i, int j, int k):
    return (i*vc.dims[1]+j)*vc.dims[2]+k

cdef inline int vc_pos_index(VolumeContainer *vc, np.float64_t *spos):
    cdef int index[3]
    cdef int i
    for i in range(3):
        index[i] = <int> ((spos[i] - vc.left_edge[i]) * vc.idds[i])
    return vc_index(vc, index[0], index[1], index[2])
