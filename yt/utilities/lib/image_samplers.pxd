"""
Definitions for image samplers




"""


import numpy as np

cimport cython
cimport numpy as np

from .partitioned_grid cimport PartitionedGrid
from .volume_container cimport VolumeContainer

DEF Nch = 4

# NOTE: We don't want to import the field_interpolator_tables here, as it
# breaks a bunch of C++ interop.  Maybe some day it won't.  So, we just forward
# declare.
cdef struct VolumeRenderAccumulator

ctypedef int calculate_extent_function(ImageSampler image,
            VolumeContainer *vc, np.int64_t rv[4]) nogil except -1

ctypedef void generate_vector_info_function(ImageSampler im,
            np.int64_t vi, np.int64_t vj,
            np.float64_t width[2],
            np.float64_t v_dir[3], np.float64_t v_pos[3]) nogil

cdef struct ImageAccumulator:
    np.float64_t rgba[Nch]
    void *supp_data

cdef class ImageSampler:
    cdef np.float64_t[:,:,:] vp_pos
    cdef np.float64_t[:,:,:] vp_dir
    cdef np.float64_t *center
    cdef np.float64_t[:,:,:] image
    cdef np.float64_t[:,:] zbuffer
    cdef np.int64_t[:,:] image_used
    cdef np.int64_t[:,:] mesh_lines
    cdef np.float64_t pdx, pdy
    cdef np.float64_t bounds[4]
    cdef np.float64_t[:,:] camera_data   # position, width, unit_vec[0,2]
    cdef int nv[2]
    cdef np.float64_t *x_vec
    cdef np.float64_t *y_vec
    cdef public object acenter, aimage, ax_vec, ay_vec
    cdef public object azbuffer
    cdef public object aimage_used
    cdef public object amesh_lines
    cdef void *supp_data
    cdef np.float64_t width[3]
    cdef public object lens_type
    cdef public str volume_method
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
    cdef object tree_containers

cdef class LightSourceRenderSampler(ImageSampler):
    cdef VolumeRenderAccumulator *vra
    cdef public object tf_obj
    cdef public object my_field_tables
