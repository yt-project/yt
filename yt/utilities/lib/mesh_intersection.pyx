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
from yt.utilities.lib.primitives cimport \
    patchSurfaceFunc, \
    patchSurfaceDerivU, \
    patchSurfaceDerivV
from vec3_ops cimport dot, subtract, cross, distance


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
                             size_t item) nogil:

    cdef Patch patch = patches[item]

    # first we compute the two planes that define the ray.
    cdef float[3] n, N1, N2
    cdef float A = dot(ray.dir, ray.dir)
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

    cdef float d1 = -dot(N1, ray.org)
    cdef float d2 = -dot(N2, ray.org)

    # the initial guess is set to zero
    cdef float u = 0.0
    cdef float v = 0.0
    cdef float[3] S
    patchSurfaceFunc(patch.v, u, v, S)
    cdef float fu = dot(N1, S) + d1
    cdef float fv = dot(N2, S) + d2
    cdef float err = fmax(fabs(fu), fabs(fv))
    
    # begin Newton interation
    cdef float tol = 1.0e-5
    cdef int iterations = 0
    cdef int max_iter = 10
    cdef float[3] Su
    cdef float[3] Sv
    cdef float J11, J12, J21, J22, det
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
    cdef float t = distance(S, ray.org)

    # only count this is it's the closest hit
    if (t < ray.tnear or t > ray.Ng[0]):
        return

    if (fabs(u) <= 1.0 and fabs(v) <= 1.0 and iterations < max_iter):

        # we have a hit, so update ray information
        ray.u = u
        ray.v = v
        ray.geomID = patch.geomID
        ray.primID = item
        ray.Ng[0] = t
        
        # sample the solution at the calculated point
        sample_hex20(patches, ray)

    return
