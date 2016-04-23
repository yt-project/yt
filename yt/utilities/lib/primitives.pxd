cimport cython 
import numpy as np
cimport numpy as np
from vec3_ops cimport dot, subtract, cross
from yt.utilities.lib.bounding_volume_hierarchy cimport Ray, BBox

cdef struct Triangle:
    np.float64_t p0[3]
    np.float64_t p1[3]
    np.float64_t p2[3]
    np.int64_t elem_id

cdef inline np.int64_t ray_triangle_intersect(const void* primitives,
                                              const np.int64_t item,
                                              Ray* ray) nogil

cdef inline void triangle_centroid(const void *primitives,
                                   const np.int64_t item,
                                   np.float64_t[3] centroid) nogil

cdef inline void triangle_bbox(const void *primitives,
                               const np.int64_t item,
                               BBox* bbox) nogil

cdef struct Patch:
    np.float64_t[8][3] v  # 8 vertices per patch
    np.int64_t elem_id
    
cdef inline np.int64_t ray_patch_intersect(const void* primitives,
                                           const np.int64_t item,
                                           Ray* ray) nogil

cdef inline void patch_centroid(const void *primitives,
                                const np.int64_t item,
                                np.float64_t[3] centroid) nogil

cdef inline void patch_bbox(const void *primitives,
                            const np.int64_t item,
                            BBox* bbox) nogil
