"""
This file contains functions used for performing ray-tracing with Embree
for 2nd-order Lagrange Elements.

Note - this file is only used for the Embree-accelerated ray-tracer.

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
from yt.utilities.lib.bounding_volume_hierarchy cimport BBox
from yt.utilities.lib.primitives cimport \
    patchSurfaceFunc, \
    patchSurfaceDerivU, \
    patchSurfaceDerivV, \
    RayHitData, \
    compute_patch_hit
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

    cdef RayHitData hd = compute_patch_hit(patch.v, ray.org, ray.dir)

    # only count this is it's the closest hit
    if (hd.t < ray.tnear or hd.t > ray.Ng[0]):
        return

    if (fabs(hd.u) <= 1.0 and fabs(hd.v) <= 1.0 and hd.converged):

        # we have a hit, so update ray information
        ray.u = hd.u
        ray.v = hd.v
        ray.geomID = patch.geomID
        ray.primID = item
        ray.Ng[0] = hd.t
        
        # sample the solution at the calculated point
        sample_hex20(patches, ray)

    return
