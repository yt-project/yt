"""
Definitions for the lens code




"""


import numpy as np

cimport cython
cimport numpy as np
from libc.math cimport (
    M_PI,
    acos,
    asin,
    atan,
    atan2,
    cos,
    exp,
    fabs,
    floor,
    log2,
    sin,
    sqrt,
)
from vec3_ops cimport L2_norm, dot, fma, subtract

from yt.utilities.lib.fp_utils cimport fclip, fmax, fmin, i64clip, iclip, imax, imin

from .image_samplers cimport (
    ImageSampler,
    calculate_extent_function,
    generate_vector_info_function,
)
from .volume_container cimport VolumeContainer


cdef extern from "platform_dep.h":
    long int lrint(double x) nogil

cdef extern from "limits.h":
    cdef int SHRT_MAX

cdef generate_vector_info_function generate_vector_info_plane_parallel
cdef generate_vector_info_function generate_vector_info_null
cdef calculate_extent_function calculate_extent_plane_parallel
cdef calculate_extent_function calculate_extent_perspective
cdef calculate_extent_function calculate_extent_null
