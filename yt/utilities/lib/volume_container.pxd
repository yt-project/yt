"""
A volume container




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np

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

cdef inline int vc_index(VolumeContainer *vc, int i, int j, int k):
    return (i*vc.dims[1]+j)*vc.dims[2]+k

cdef inline int vc_pos_index(VolumeContainer *vc, np.float64_t *spos):
    cdef int index[3]
    cdef int i
    for i in range(3):
        index[i] = <int> ((spos[i] - vc.left_edge[i]) * vc.idds[i])
    return vc_index(vc, index[0], index[1], index[2])

