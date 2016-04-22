cimport cython 
import numpy as np
cimport numpy as np
from vec3_ops cimport dot, subtract, cross
from yt.utilities.lib.bounding_volume_hierarchy cimport Ray, BBox

cdef np.float64_t DETERMINANT_EPS = 1.0e-10

cdef extern from "platform_dep.h" nogil:
    double fmax(double x, double y)
    double fmin(double x, double y)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline np.int64_t ray_triangle_intersect(const void* primitives,
                                              const np.int64_t item,
                                              Ray* ray) nogil:
# https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm

    cdef Triangle tri = (<Triangle*> primitives)[item]

    # edge vectors
    cdef np.float64_t e1[3]
    cdef np.float64_t e2[3]
    subtract(tri.p1, tri.p0, e1)
    subtract(tri.p2, tri.p0, e2)

    cdef np.float64_t P[3]
    cross(ray.direction, e2, P)

    cdef np.float64_t det, inv_det
    det = dot(e1, P)
    if(det > -DETERMINANT_EPS and det < DETERMINANT_EPS): 
        return False
    inv_det = 1.0 / det

    cdef np.float64_t T[3]
    subtract(ray.origin, tri.p0, T)

    cdef np.float64_t u = dot(T, P) * inv_det
    if(u < 0.0 or u > 1.0):
        return False

    cdef np.float64_t Q[3]
    cross(T, e1, Q)

    cdef np.float64_t v = dot(ray.direction, Q) * inv_det
    if(v < 0.0 or u + v  > 1.0):
        return False

    cdef np.float64_t t = dot(e2, Q) * inv_det

    if(t > DETERMINANT_EPS and t < ray.t_far):
        ray.t_far = t
        ray.elem_id = tri.elem_id
        return True

    return False


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void triangle_centroid(const void *primitives,
                                   const np.int64_t item,
                                   np.float64_t[3] centroid) nogil:

    cdef Triangle tri = (<Triangle*> primitives)[item]
    cdef np.int64_t i
    for i in range(3):
        centroid[i] = (tri.p0[i] + tri.p1[i] + tri.p2[i]) / 3.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void triangle_bbox(const void *primitives,
                               const np.int64_t item,
                               BBox* bbox) nogil:

    cdef Triangle tri = (<Triangle*> primitives)[item]
    cdef np.int64_t i
    for i in range(3):
        bbox.left_edge[i] = fmin(fmin(tri.p0[i], tri.p1[i]), tri.p2[i])
        bbox.right_edge[i] = fmax(fmax(tri.p0[i], tri.p1[i]), tri.p2[i])
