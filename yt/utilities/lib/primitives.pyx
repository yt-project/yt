cimport cython
import numpy as np
cimport numpy as np
cimport cython.floating
from libc.math cimport fabs

from yt.utilities.lib.vec3_ops cimport dot, subtract, cross, distance, L2_norm

cdef np.float64_t DETERMINANT_EPS = 1.0e-10
cdef np.float64_t INF = np.inf

cdef extern from "platform_dep.h" nogil:
    double fmax(double x, double y)
    double fmin(double x, double y)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.int64_t ray_bbox_intersect(Ray* ray, const BBox bbox) nogil:
    '''
    
    This returns an integer flag that indicates whether a ray and a bounding
    box intersect. It does not modify either either the ray or the box.
    
    '''

    # https://tavianator.com/fast-branchless-raybounding-box-intersections/

    cdef np.float64_t tmin = -INF
    cdef np.float64_t tmax =  INF
 
    cdef np.int64_t i
    cdef np.float64_t t1, t2
    for i in range(3):
        t1 = (bbox.left_edge[i]  - ray.origin[i])*ray.inv_dir[i]
        t2 = (bbox.right_edge[i] - ray.origin[i])*ray.inv_dir[i] 
        tmin = fmax(tmin, fmin(t1, t2))
        tmax = fmin(tmax, fmax(t1, t2))
 
    return tmax >= fmax(tmin, 0.0)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.int64_t ray_triangle_intersect(const void* primitives,
                                       const np.int64_t item,
                                       Ray* ray) nogil:
    '''
    
    This returns an integer flag that indicates whether a triangle is the
    closest hit for the ray so far. If it is, the ray is updated to store the
    current triangle index and the distance to the first hit. The triangle used
    is the one indexed by "item" in the array of primitives.
    
    
    '''

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
cdef void triangle_centroid(const void *primitives,
                            const np.int64_t item,
                            np.float64_t[3] centroid) nogil:
    '''
    
    This computes the centroid of the input triangle. The triangle used
    is the one indexed by "item" in the array of primitives. The result
    will be stored in the numpy array passed in as "centroid".
    
    '''
        
    cdef Triangle tri = (<Triangle*> primitives)[item]
    cdef np.int64_t i
    for i in range(3):
        centroid[i] = (tri.p0[i] + tri.p1[i] + tri.p2[i]) / 3.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void triangle_bbox(const void *primitives,
                        const np.int64_t item,
                        BBox* bbox) nogil:
    '''
    
    This computes the bounding box of the input triangle. The triangle used
    is the one indexed by "item" in the array of primitives. The result
    will be stored in the input BBox.
    
    '''
    
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
    '''
    
    This function is a parametric representation of the surface of a bi-quadratic
    patch. The inputs are the eight nodes that define a face of a 20-node hex element,
    and two parameters u and v that vary from -1 to 1 and tell you where you are on
    the surface of the patch. The output is the array 'S' that stores the physical
    (x, y, z) position of the corresponding point on the patch. This function is needed
    to compute the intersection of rays and bi-quadratic patches.
    
    '''
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
    '''
    
    This function computes the derivative of the S(u, v) function w.r.t u. 
    
    '''
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
    '''
    
    This function computes the derivative of the S(u, v) function w.r.t v.
    
    '''
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
cdef RayHitData compute_patch_hit(cython.floating[8][3] verts,
                                  cython.floating[3] ray_origin,
                                  cython.floating[3] ray_direction) nogil:
    """
    
    This function iteratively computes whether the bi-quadratic patch defined by the
    eight input nodes intersects with the given ray. Either way, information about
    the potential hit is stored in the returned RayHitData.
    
    """
    # first we compute the two planes that define the ray.
    cdef cython.floating[3] n, N1, N2
    cdef cython.floating A = dot(ray_direction, ray_direction)
    for i in range(3):
        n[i] = ray_direction[i] / A

    if ((fabs(n[0]) > fabs(n[1])) and (fabs(n[0]) > fabs(n[2]))):
        N1[0] = n[1]
        N1[1] =-n[0]
        N1[2] = 0.0
    else:
        N1[0] = 0.0
        N1[1] = n[2]
        N1[2] =-n[1]
    cross(N1, n, N2)

    cdef cython.floating d1 = -dot(N1, ray_origin)
    cdef cython.floating d2 = -dot(N2, ray_origin)

    # the initial guess is set to zero
    cdef cython.floating u = 0.0
    cdef cython.floating v = 0.0
    cdef cython.floating[3] S
    patchSurfaceFunc(verts, u, v, S)
    cdef cython.floating fu = dot(N1, S) + d1
    cdef cython.floating fv = dot(N2, S) + d2
    cdef cython.floating err = fmax(fabs(fu), fabs(fv))
    
    # begin Newton interation
    cdef cython.floating tol = 1.0e-5
    cdef int iterations = 0
    cdef int max_iter = 10
    cdef cython.floating[3] Su
    cdef cython.floating[3] Sv
    cdef cython.floating J11, J12, J21, J22, det
    while ((err > tol) and (iterations < max_iter)):
        # compute the Jacobian
        patchSurfaceDerivU(verts, u, v, Su)
        patchSurfaceDerivV(verts, u, v, Sv)
        J11 = dot(N1, Su)
        J12 = dot(N1, Sv)
        J21 = dot(N2, Su)
        J22 = dot(N2, Sv)
        det = (J11*J22 - J12*J21)
        
        # update the u, v values
        u -= ( J22*fu - J12*fv) / det
        v -= (-J21*fu + J11*fv) / det
        
        patchSurfaceFunc(verts, u, v, S)
        fu = dot(N1, S) + d1
        fv = dot(N2, S) + d2

        err = fmax(fabs(fu), fabs(fv))
        iterations += 1

    # t is the distance along the ray to this hit
    cdef cython.floating t = distance(S, ray_origin) / L2_norm(ray_direction)
    
    # return hit data
    cdef RayHitData hd
    hd.u = u
    hd.v = v
    hd.t = t
    hd.converged = (iterations < max_iter)
    return hd


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.int64_t ray_patch_intersect(const void* primitives,
                                    const np.int64_t item,
                                    Ray* ray) nogil:
    '''

    This returns an integer flag that indicates whether the given patch is the
    closest hit for the ray so far. If it is, the ray is updated to store the
    current primitive index and the distance to the first hit. The patch used
    is the one indexed by "item" in the array of primitives.


    '''
    cdef Patch patch = (<Patch*> primitives)[item]

    cdef RayHitData hd = compute_patch_hit(patch.v, ray.origin, ray.direction)
    
    # only count this is it's the closest hit
    if (hd.t < ray.t_near or hd.t > ray.t_far):
        return False

    if (fabs(hd.u) <= 1.0 and fabs(hd.v) <= 1.0 and hd.converged):
        # we have a hit, so update ray information
        ray.t_far = hd.t
        ray.elem_id = patch.elem_id
        return True

    return False
        

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void patch_centroid(const void *primitives,
                         const np.int64_t item,
                         np.float64_t[3] centroid) nogil:
    '''

    This computes the centroid of the input patch. The patch used
    is the one indexed by "item" in the array of primitives. The result
    will be stored in the numpy array passed in as "centroid".

    '''

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
cdef void patch_bbox(const void *primitives,
                    const np.int64_t item,
                     BBox* bbox) nogil:

    '''

    This computes the bounding box of the input patch. The patch used
    is the one indexed by "item" in the array of primitives. The result
    will be stored in the input BBox.

    '''

    cdef np.int64_t i, j
    cdef Patch patch = (<Patch*> primitives)[item]
    
    for j in range(3):
        bbox.left_edge[j] = patch.v[0][j]
        bbox.right_edge[j] = patch.v[0][j]

    for i in range(1, 8):
        for j in range(3):
            bbox.left_edge[j] = fmin(bbox.left_edge[j], patch.v[i][j])
            bbox.right_edge[j] = fmax(bbox.right_edge[j], patch.v[i][j])
