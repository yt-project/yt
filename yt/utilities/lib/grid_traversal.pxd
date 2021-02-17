"""
Definitions for the traversal code




"""


import numpy as np

cimport cython
cimport numpy as np

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

#-----------------------------------------------------------------------------
# walk_volume(VolumeContainer *vc,  np.float64_t v_pos[3], np.float64_t v_dir[3], sampler_function *sample,
#             void *data, np.float64_t *return_t = NULL, np.float64_t max_t = 1.0)
#    vc        VolumeContainer*  : Pointer to the volume container to be traversed.
#    v_pos     np.float64_t[3]   : The x,y,z coordinates of the ray's origin.
#    v_dir     np.float64_t[3]   : The x,y,z coordinates of the ray's direction.
#    sample    sampler_function* : Pointer to the sampler function to be used.
#    return_t  np.float64_t*     : Pointer to the final value of t that is still inside the volume container. Defaulted to NULL.
#    max_t     np.float64_t      : The maximum value of t that the ray is allowed to travel. Defaulted to 1.0 (no restriction).
#
#    Note: 't' is not time here. Rather, it is a factor representing the difference between the initial point 'v_pos'
#             and the end point, which we might call v_end. It is scaled such that v_pos + v * t = v_pos at t = 0.0, and
#             v_end at t = 1.0. Therefore, if max_t is set to 1.0, there is no restriction on t.
#
# Written by the yt Development Team.
# Encapsulates the Amanatides & Woo "Fast Traversal Voxel Algorithm" to walk over a volume container 'vc'
# The function occurs in two phases, initialization and traversal.
# See: https://www.researchgate.net/publication/2611491_A_Fast_Voxel_Traversal_Algorithm_for_Ray_Tracing
# Returns: The number of voxels hit during the traversal phase. If the traversal phase is not reached, returns 0.
#-----------------------------------------------------------------------------
cdef int walk_volume(VolumeContainer *vc,
                     np.float64_t v_pos[3],
                     np.float64_t v_dir[3],
                     sampler_function *sampler,
                     void *data,
                     np.float64_t *return_t = *,
                     np.float64_t max_t = *) nogil
