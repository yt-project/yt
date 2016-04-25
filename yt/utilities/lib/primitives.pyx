cimport cython 
import numpy as np
cimport numpy as np
cimport cython.floating
from libc.math cimport fabs

from yt.utilities.lib.vec3_ops cimport dot, subtract, cross, distance, L2_norm
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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void patchSurfaceFunc(const cython.floating[8][3] verts,
                           const cython.floating u, 
                           const cython.floating v,
                           cython.floating[3] S) nogil:

  cdef int i
  for i in range(3):
      S[i] = 0.25*(1.0 - u)*(1.0 - v)*(-u - v - 1)*verts[0][i] + \
             0.25*(1.0 + u)*(1.0 - v)*( u - v - 1)*verts[1][i] + \
             0.25*(1.0 + u)*(1.0 + v)*( u + v - 1)*verts[2][i] + \
             0.25*(1.0 - u)*(1.0 + v)*(-u + v - 1)*verts[3][i] + \
             0.5*(1 - u)*(1 - v*v)*verts[4][i] + \
             0.5*(1 - u*u)*(1 - v)*verts[5][i] + \
             0.5*(1 + u)*(1 - v*v)*verts[6][i] + \
             0.5*(1 - u*u)*(1 + v)*verts[7][i]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void patchSurfaceDerivU(const cython.floating[8][3] verts,
                             const cython.floating u, 
                             const cython.floating v,
                             cython.floating[3] Su) nogil: 
  cdef int i
  for i in range(3):
      Su[i] = (-0.25*(v - 1.0)*(u + v + 1) - 0.25*(u - 1.0)*(v - 1.0))*verts[0][i] + \
              (-0.25*(v - 1.0)*(u - v - 1) - 0.25*(u + 1.0)*(v - 1.0))*verts[1][i] + \
              ( 0.25*(v + 1.0)*(u + v - 1) + 0.25*(u + 1.0)*(v + 1.0))*verts[2][i] + \
              ( 0.25*(v + 1.0)*(u - v + 1) + 0.25*(u - 1.0)*(v + 1.0))*verts[3][i] + \
              0.5*(v*v - 1.0)*verts[4][i] + u*(v - 1.0)*verts[5][i] - \
              0.5*(v*v - 1.0)*verts[6][i] - u*(v + 1.0)*verts[7][i]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void patchSurfaceDerivV(const cython.floating[8][3] verts,
                             const cython.floating u, 
                             const cython.floating v,
                             cython.floating[3] Sv) nogil:
    cdef int i 
    for i in range(3):
        Sv[i] = (-0.25*(u - 1.0)*(u + v + 1) - 0.25*(u - 1.0)*(v - 1.0))*verts[0][i] + \
                (-0.25*(u + 1.0)*(u - v - 1) + 0.25*(u + 1.0)*(v - 1.0))*verts[1][i] + \
                ( 0.25*(u + 1.0)*(u + v - 1) + 0.25*(u + 1.0)*(v + 1.0))*verts[2][i] + \
                ( 0.25*(u - 1.0)*(u - v + 1) - 0.25*(u - 1.0)*(v + 1.0))*verts[3][i] + \
                0.5*(u*u - 1.0)*verts[5][i] + v*(u - 1.0)*verts[4][i] - \
                0.5*(u*u - 1.0)*verts[7][i] - v*(u + 1.0)*verts[6][i]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline np.int64_t ray_patch_intersect(const void* primitives,
                                           const np.int64_t item,
                                           Ray* ray) nogil:

    cdef Patch patch = (<Patch*> primitives)[item]

    # first we compute the two planes that define the ray.
    cdef np.float64_t[3] n, N1, N2
    cdef np.float64_t A = dot(ray.direction, ray.direction)
    for i in range(3):
        n[i] = ray.direction[i] / A

    if ((fabs(n[0]) > fabs(n[1])) and (fabs(n[0]) > fabs(n[2]))):
        N1[0] = n[1]
        N1[1] =-n[0]
        N1[2] = 0.0
    else:
        N1[0] = 0.0
        N1[1] = n[2]
        N1[2] =-n[1]
    cross(N1, n, N2)

    cdef np.float64_t d1 = -dot(N1, ray.origin)
    cdef np.float64_t d2 = -dot(N2, ray.origin)

    # the initial guess is set to zero
    cdef np.float64_t u = 0.0
    cdef np.float64_t v = 0.0
    cdef np.float64_t[3] S
    patchSurfaceFunc(patch.v, u, v, S)
    cdef np.float64_t fu = dot(N1, S) + d1
    cdef np.float64_t fv = dot(N2, S) + d2
    cdef np.float64_t err = fmax(fabs(fu), fabs(fv))
    
    # begin Newton interation
    cdef np.float64_t tol = 1.0e-5
    cdef int iterations = 0
    cdef int max_iter = 10
    cdef np.float64_t[3] Su
    cdef np.float64_t[3] Sv
    cdef np.float64_t J11, J12, J21, J22, det
    while ((err > tol) and (iterations < max_iter)):
        # compute the Jacobian
        patchSurfaceDerivU(patch.v, u, v, Su)
        patchSurfaceDerivV(patch.v, u, v, Sv)
        J11 = dot(N1, Su)
        J12 = dot(N1, Sv)
        J21 = dot(N2, Su)
        J22 = dot(N2, Sv)
        det = (J11*J22 - J12*J21)
        
        # update the u, v values
        u -= ( J22*fu - J12*fv) / det
        v -= (-J21*fu + J11*fv) / det
        
        patchSurfaceFunc(patch.v, u, v, S)
        fu = dot(N1, S) + d1
        fv = dot(N2, S) + d2

        err = fmax(fabs(fu), fabs(fv))
        iterations += 1

    # t is the distance along the ray to this hit
    cdef np.float64_t t = distance(S, ray.origin) / L2_norm(ray.direction)

    # only count this is it's the closest hit
    if (t < ray.t_near or t > ray.t_far):
        return False

    if (fabs(u) <= 1.0 and fabs(v) <= 1.0 and iterations < max_iter):

        # we have a hit, so update ray information
        ray.t_far = t
        ray.elem_id = patch.elem_id
        return True

    return False
        

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void patch_centroid(const void *primitives,
                                const np.int64_t item,
                                np.float64_t[3] centroid) nogil:

    cdef np.int64_t i, j
    cdef Patch patch = (<Patch*> primitives)[item]

    for j in range(3):
        centroid[j] = 0.0

    for i in range(8):
        for j in range(3):
            centroid[j] += patch.v[i][j]

    for j in range(3):
        centroid[j] /= 8.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void patch_bbox(const void *primitives,
                            const np.int64_t item,
                            BBox* bbox) nogil:

    cdef np.int64_t i, j
    cdef Patch patch = (<Patch*> primitives)[item]
    
    for j in range(3):
        bbox.left_edge[j] = patch.v[0][j]
        bbox.right_edge[j] = patch.v[0][j]

    for i in range(1, 8):
        for j in range(3):
            bbox.left_edge[j] = fmin(bbox.left_edge[j], patch.v[i][j])
            bbox.right_edge[j] = fmax(bbox.right_edge[j], patch.v[i][j])
