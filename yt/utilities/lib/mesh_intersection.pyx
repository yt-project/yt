"""
This file contains functions used for performing ray-tracing with 2nd-order Lagrange
Elements.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr
cimport pyembree.rtcore_geometry as rtcg
from pyembree.rtcore cimport Vec3f
cimport numpy as np
cimport cython
from libc.math cimport fabs, fmin, fmax, sqrt
from yt.utilities.lib.mesh_samplers cimport sample_hex20

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void patchSurfaceFunc(const Patch& patch, 
                           const float u, 
                           const float v,
                           float[3] S) nogil:

  cdef int i
  for i in range(3):
      S[i] = 0.25*(1.0 - u)*(1.0 - v)*(-u - v - 1)*patch.v[0][i] + \
             0.25*(1.0 + u)*(1.0 - v)*( u - v - 1)*patch.v[1][i] + \
             0.25*(1.0 + u)*(1.0 + v)*( u + v - 1)*patch.v[2][i] + \
             0.25*(1.0 - u)*(1.0 + v)*(-u + v - 1)*patch.v[3][i] + \
             0.5*(1 - u)*(1 - v*v)*patch.v[4][i] + \
             0.5*(1 - u*u)*(1 - v)*patch.v[5][i] + \
             0.5*(1 + u)*(1 - v*v)*patch.v[6][i] + \
             0.5*(1 - u*u)*(1 + v)*patch.v[7][i]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void patchSurfaceDerivU(const Patch& patch, 
                             const float u, 
                             const float v,
                             float[3] Su) nogil: 
  cdef int i
  for i in range(3):
      Su[i] = (-0.25*(v - 1.0)*(u + v + 1) - 0.25*(u - 1.0)*(v - 1.0))*patch.v[0][i] + \
              (-0.25*(v - 1.0)*(u - v - 1) - 0.25*(u + 1.0)*(v - 1.0))*patch.v[1][i] + \
              ( 0.25*(v + 1.0)*(u + v - 1) + 0.25*(u + 1.0)*(v + 1.0))*patch.v[2][i] + \
              ( 0.25*(v + 1.0)*(u - v + 1) + 0.25*(u - 1.0)*(v + 1.0))*patch.v[3][i] + \
              0.5*(v*v - 1.0)*patch.v[4][i] + u*(v - 1.0)*patch.v[5][i] - \
              0.5*(v*v - 1.0)*patch.v[6][i] - u*(v + 1.0)*patch.v[7][i]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void patchSurfaceDerivV(const Patch& patch, 
                             const float u, 
                             const float v,
                             float[3] Sv) nogil:
    cdef int i 
    for i in range(3):
        Sv[i] = (-0.25*(u - 1.0)*(u + v + 1) - 0.25*(u - 1.0)*(v - 1.0))*patch.v[0][i] + \
                (-0.25*(u + 1.0)*(u - v - 1) + 0.25*(u + 1.0)*(v - 1.0))*patch.v[1][i] + \
                ( 0.25*(u + 1.0)*(u + v - 1) + 0.25*(u + 1.0)*(v + 1.0))*patch.v[2][i] + \
                ( 0.25*(u - 1.0)*(u - v + 1) - 0.25*(u - 1.0)*(v + 1.0))*patch.v[3][i] + \
                0.5*(u*u - 1.0)*patch.v[5][i] + v*(u - 1.0)*patch.v[4][i] - \
                0.5*(u*u - 1.0)*patch.v[7][i] - v*(u + 1.0)*patch.v[6][i]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef float dot(const float* a, 
               const float* b,
               size_t N) nogil:
    cdef int i
    cdef float rv = 0.0
    for i in range(N):
        rv += a[i]*b[i]
    return rv


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void cross(const float* a, 
                const float* b,
                float* c) nogil:
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void patchBoundsFunc(Patch* patches, 
                          size_t item, 
                          rtcg.RTCBounds* bounds_o) nogil:

    cdef Patch patch = patches[item]
    
    cdef float lo_x = 1.0e300
    cdef float lo_y = 1.0e300
    cdef float lo_z = 1.0e300

    cdef float hi_x = -1.0e300
    cdef float hi_y = -1.0e300
    cdef float hi_z = -1.0e300

    cdef int i
    for i in range(8):
        lo_x = fmin(lo_x, patch.v[i][0])
        lo_y = fmin(lo_y, patch.v[i][1])
        lo_z = fmin(lo_z, patch.v[i][2])
        hi_x = fmax(hi_x, patch.v[i][0])
        hi_y = fmax(hi_y, patch.v[i][1])
        hi_z = fmax(hi_z, patch.v[i][2])

    bounds_o.lower_x = lo_x
    bounds_o.lower_y = lo_y
    bounds_o.lower_z = lo_z
    bounds_o.upper_x = hi_x
    bounds_o.upper_y = hi_y
    bounds_o.upper_z = hi_z


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void patchIntersectFunc(Patch* patches,
                             rtcr.RTCRay& ray,
                             size_t item):

    cdef Patch patch = patches[item]

    # return if the ray does not intersect the bounding sphere
    cdef float[3] vec
    for i in range(3):
        vec[i] = ray.org[i] - patch.center[i]
        
    cdef float A = dot(ray.dir, ray.dir, 3)
    cdef float B = 2.0*dot(vec, ray.dir, 3)
    cdef float C = dot(vec, vec, 3) - patch.radius**2
    cdef float D = B*B - 4.0*A*C
    if (D < 0.0):
        return

    # otherwise, iterate to the get the true hit position
    # first we compute the the two planes that define the ray.
    cdef float[3] n, N1, N2
    for i in range(3):
        n[i] = ray.dir[i] / A

    if ((fabs(n[0]) > fabs(n[1])) and (fabs(n[0]) > fabs(n[2]))):
        N1[0] = n[1]
        N1[1] =-n[0]
        N1[2] = 0.0
    else:
        N1[0] = 0.0
        N1[1] = n[2]
        N1[2] =-n[1]
    cross(N1, n, N2)

    cdef float d1 = -dot(N1, ray.org, 3)
    cdef float d2 = -dot(N2, ray.org, 3)

    # the initial guess is set to zero
    cdef float u = 0.0
    cdef float v = 0.0
    cdef float[3] S
    patchSurfaceFunc(patch, u, v, S)
    cdef float fu = dot(N1, S, 3) + d1
    cdef float fv = dot(N2, S, 3) + d2
    cdef float err = fmax(fabs(fu), fabs(fv))
    
    # begin Newton interation
    cdef float tol = 1.0e-6
    cdef int iterations = 0
    cdef int max_iter = 10
    cdef float[3] Su
    cdef float[3] Sv
    cdef float J11, J12, J21, J22, det
    while ((err > tol) and (iterations < max_iter)):
        # compute the Jacobian
        patchSurfaceDerivU(patch, u, v, Su)
        patchSurfaceDerivV(patch, u, v, Sv)
        J11 = dot(N1, Su, 3)
        J12 = dot(N1, Sv, 3)
        J21 = dot(N2, Su, 3)
        J22 = dot(N2, Sv, 3)
        det = (J11*J22 - J12*J21)
        
        # update the u, v values
        u -= ( J22*fu - J12*fv) / det
        v -= (-J21*fu + J11*fv) / det
        
        patchSurfaceFunc(patch, u, v, S)
        fu = dot(N1, S, 3) + d1
        fv = dot(N2, S, 3) + d2

        err = fmax(fabs(fu), fabs(fv))
        iterations += 1

    # only count this as a hit if it's within the current ray endpoints
    cdef float t = 0.0
    for i in range(3):
        t += (S[i] - ray.org[i])**2
    t = sqrt(t)
    if (t < ray.tnear or t > ray.tfar):
        return

    patchSurfaceDerivU(patch, u, v, Su)
    patchSurfaceDerivV(patch, u, v, Sv)
    
    if (fabs(u) <= 1.0 and fabs(v) <= 1.0 and iterations < max_iter):
        # we have a hit, so update ray information
        ray.u = u
        ray.v = v
        ray.tfar = t
        ray.geomID = patch.geomID
        ray.primID = item
        cross(Su, Sv, ray.Ng)

        # sample the solution at the calculated point
        sample_hex20(patches, ray)
        ray.time = (item % 6)
    return
