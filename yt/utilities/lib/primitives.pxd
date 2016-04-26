cimport cython 
cimport cython.floating
import numpy as np
cimport numpy as np
from vec3_ops cimport dot, subtract, cross
from yt.utilities.lib.bounding_volume_hierarchy cimport Ray, BBox

cdef struct RayHitData:
    np.float64_t u
    np.float64_t v
    np.float64_t t
    np.int64_t converged

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

cdef void patchSurfaceFunc(const cython.floating[8][3] verts,
                           const cython.floating u,
                           const cython.floating v,
                           cython.floating[3] S) nogil

cdef void patchSurfaceDerivU(const cython.floating[8][3] verts,
                             const cython.floating u,
                             const cython.floating v,
                             cython.floating[3] Su) nogil

cdef void patchSurfaceDerivV(const cython.floating[8][3] verts,
                             const cython.floating u,
                             const cython.floating v,
                             cython.floating[3] Sv) nogil

cdef RayHitData compute_patch_hit(cython.floating[8][3] verts,
                                  cython.floating[3] ray_origin,
                                  cython.floating[3] ray_direction) nogil
    
cdef inline np.int64_t ray_patch_intersect(const void* primitives,
                                           const np.int64_t item,
                                           Ray* ray) nogil

cdef inline void patch_centroid(const void *primitives,
                                const np.int64_t item,
                                np.float64_t[3] centroid) nogil

cdef inline void patch_bbox(const void *primitives,
                            const np.int64_t item,
                            BBox* bbox) nogil
