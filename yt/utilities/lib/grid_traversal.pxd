"""
Definitions for the traversal code




"""


import numpy as np
cimport numpy as np
cimport cython
from .image_samplers cimport ImageSampler
from .volume_container cimport VolumeContainer, vc_index, vc_pos_index

ctypedef void sampler_function(
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
                     sampler_function *sampler,
                     void *data,
                     np.float64_t *return_t = *,
                     np.float64_t max_t = *) nogil

