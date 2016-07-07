"""
Definitions for image samplers




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
cimport kdtree_utils
from .volume_container cimport VolumeContainer
from .lenses cimport \
    calculate_extent_function, \
    generate_vector_info_function
from .partitioned_grid cimport PartitionedGrid

from field_interpolation_tables cimport \
    FieldInterpolationTable

DEF Nch = 4

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

cdef struct VolumeRenderAccumulator:
    int n_fits
    int n_samples
    FieldInterpolationTable *fits
    int field_table_ids[6]
    np.float64_t star_coeff
    np.float64_t star_er
    np.float64_t star_sigma_num
    kdtree_utils.kdtree *star_list
    np.float64_t *light_dir
    np.float64_t *light_rgba
    int grey_opacity

cdef struct ImageAccumulator:
    np.float64_t rgba[Nch]
    void *supp_data

cdef class ImageSampler:
    cdef ImageContainer *image
    cdef public object acenter, aimage, ax_vec, ay_vec
    cdef public object azbuffer
    cdef public object aimage_used
    cdef public object amesh_lines
    cdef void *supp_data
    cdef np.float64_t width[3]
    cdef public object lens_type
    cdef calculate_extent_function *extent_function
    cdef generate_vector_info_function *vector_function
    cdef void setup(self, PartitionedGrid pg)
    @staticmethod
    cdef void sample(VolumeContainer *vc,
                np.float64_t v_pos[3],
                np.float64_t v_dir[3],
                np.float64_t enter_t,
                np.float64_t exit_t,
                int index[3],
                void *data) nogil

cdef class ProjectionSampler(ImageSampler):
    pass

cdef class InterpolatedProjectionSampler(ImageSampler):
    cdef VolumeRenderAccumulator *vra
    cdef public object tf_obj
    cdef public object my_field_tables

cdef class VolumeRenderSampler(ImageSampler):
    cdef VolumeRenderAccumulator *vra
    cdef public object tf_obj
    cdef public object my_field_tables
    cdef kdtree_utils.kdtree **trees
    cdef object tree_containers

cdef class LightSourceRenderSampler(ImageSampler):
    cdef VolumeRenderAccumulator *vra
    cdef public object tf_obj
    cdef public object my_field_tables
