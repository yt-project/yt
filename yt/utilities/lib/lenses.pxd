"""
Definitions for the lens code




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython
from .volume_container cimport VolumeContainer
from vec3_ops cimport dot, subtract, L2_norm, fma
from libc.math cimport exp, floor, log2, \
    fabs, atan, atan2, asin, cos, sin, sqrt, acos, M_PI
from yt.utilities.lib.fp_utils cimport imax, fmax, imin, fmin, iclip, fclip, i64clip

cdef extern from "platform_dep.h":
    long int lrint(double x) nogil

cdef extern from "limits.h":
    cdef int SHRT_MAX

cdef struct ImageContainer:
    np.float64_t[:,:,:] vp_pos
    np.float64_t[:,:,:] vp_dir
    np.float64_t *center
    np.float64_t[:,:,:] image
    np.float64_t[:,:] zbuffer
    np.int64_t[:,:] image_used
    np.int64_t[:,:] mesh_lines
    np.float64_t pdx, pdy
    np.float64_t bounds[4]
    np.float64_t[:,:] camera_data   # position, width, unit_vec[0,2]
    int nv[2]
    np.float64_t *x_vec
    np.float64_t *y_vec


ctypedef int calculate_extent_function(ImageContainer *image,
            VolumeContainer *vc, np.int64_t rv[4]) nogil except -1

ctypedef void generate_vector_info_function(ImageContainer *im,
            np.int64_t vi, np.int64_t vj,
            np.float64_t width[2],
            np.float64_t v_dir[3], np.float64_t v_pos[3]) nogil

cdef generate_vector_info_function generate_vector_info_plane_parallel
cdef generate_vector_info_function generate_vector_info_null
cdef calculate_extent_function calculate_extent_plane_parallel
cdef calculate_extent_function calculate_extent_perspective
cdef calculate_extent_function calculate_extent_null
