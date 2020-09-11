"""
Definitions for the lens code




"""


import numpy as np
cimport numpy as np
cimport cython
from .volume_container cimport VolumeContainer
from vec3_ops cimport dot, subtract, L2_norm, fma
from libc.math cimport exp, floor, log2, \
    fabs, atan, atan2, asin, cos, sin, sqrt, acos, M_PI
from yt.utilities.lib.fp_utils cimport imax, fmax, imin, fmin, iclip, fclip, i64clip
from .image_samplers cimport \
    ImageSampler, \
    calculate_extent_function, \
    generate_vector_info_function

cdef extern from "platform_dep.h":
    long int lrint(double x) nogil

cdef extern from "limits.h":
    cdef int SHRT_MAX

cdef generate_vector_info_function generate_vector_info_plane_parallel
cdef generate_vector_info_function generate_vector_info_null
cdef calculate_extent_function calculate_extent_plane_parallel
cdef calculate_extent_function calculate_extent_perspective
cdef calculate_extent_function calculate_extent_null
