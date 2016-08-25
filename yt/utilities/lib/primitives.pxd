cimport cython 
cimport cython.floating
import numpy as np
cimport numpy as np
from vec3_ops cimport dot, subtract, cross

cdef struct Ray:
    np.float64_t origin[3]
    np.float64_t direction[3]
    np.float64_t inv_dir[3]
    np.float64_t data_val
    np.float64_t t_near
    np.float64_t t_far
    np.int64_t elem_id
    np.int64_t near_boundary

cdef struct BBox:
    np.float64_t left_edge[3]
    np.float64_t right_edge[3]

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

cdef np.int64_t ray_bbox_intersect(Ray* ray, const BBox bbox) nogil

cdef np.int64_t ray_triangle_intersect(const void* primitives,
                                       const np.int64_t item,
                                       Ray* ray) nogil

cdef void triangle_centroid(const void *primitives,
                            const np.int64_t item,
                            np.float64_t[3] centroid) nogil

cdef void triangle_bbox(const void *primitives,
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
    
cdef np.int64_t ray_patch_intersect(const void* primitives,
                                    const np.int64_t item,
                                    Ray* ray) nogil

cdef void patch_centroid(const void *primitives,
                         const np.int64_t item,
                         np.float64_t[3] centroid) nogil

cdef void patch_bbox(const void *primitives,
                     const np.int64_t item,
                     BBox* bbox) nogil
