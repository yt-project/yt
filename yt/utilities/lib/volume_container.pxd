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
cimport cython
cimport cython.view

cdef struct VolumeContainer:
    int n_fields
    np.float64_t[::cython.view.indirect,:,:,::1]  data
    # The mask has dimensions one fewer in each direction than data
    np.uint8_t[:,:,:] mask
    np.float64_t left_edge[3]
    np.float64_t right_edge[3]
    np.float64_t dds[3]
    np.float64_t idds[3]
    int dims[3]
