# distutils: include_dirs = EMBREE_INC_DIR
# distutils: library_dirs = EMBREE_LIB_DIR
# distutils: libraries = EMBREE_LIBS
# distutils: language = c++
"""
This file contains functions used for performing ray-tracing with Embree
for 2nd-order Lagrange Elements.

Note - this file is only used for the Embree-accelerated ray-tracer.

"""


cimport cython
cimport numpy as np
cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_geometry as rtcg
cimport pyembree.rtcore_ray as rtcr
from libc.math cimport fabs, fmax, fmin, sqrt
from pyembree.rtcore cimport Vec3f

from yt.utilities.lib.bounding_volume_hierarchy cimport BBox
from yt.utilities.lib.primitives cimport (
    RayHitData,
    compute_patch_hit,
    compute_tet_patch_hit,
    patchSurfaceDerivU,
    patchSurfaceDerivV,
    patchSurfaceFunc,
    tet_patchSurfaceDerivU,
    tet_patchSurfaceDerivV,
    tet_patchSurfaceFunc,
)
from yt.utilities.lib.vec3_ops cimport cross, distance, dot, subtract

from .mesh_samplers cimport sample_hex20, sample_tet10


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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void tet_patchBoundsFunc(Tet_Patch* tet_patches,
                          size_t item,
                          rtcg.RTCBounds* bounds_o) nogil:

    cdef Tet_Patch tet_patch = tet_patches[item]

    cdef float lo_x = 1.0e300
    cdef float lo_y = 1.0e300
    cdef float lo_z = 1.0e300

    cdef float hi_x = -1.0e300
    cdef float hi_y = -1.0e300
    cdef float hi_z = -1.0e300

    cdef int i
    for i in range(6):
        lo_x = fmin(lo_x, tet_patch.v[i][0])
        lo_y = fmin(lo_y, tet_patch.v[i][1])
        lo_z = fmin(lo_z, tet_patch.v[i][2])
        hi_x = fmax(hi_x, tet_patch.v[i][0])
        hi_y = fmax(hi_y, tet_patch.v[i][1])
        hi_z = fmax(hi_z, tet_patch.v[i][2])

    bounds_o.lower_x = lo_x
    bounds_o.lower_y = lo_y
    bounds_o.lower_z = lo_z
    bounds_o.upper_x = hi_x
    bounds_o.upper_y = hi_y
    bounds_o.upper_z = hi_z


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void tet_patchIntersectFunc(Tet_Patch* tet_patches,
                             rtcr.RTCRay& ray,
                             size_t item) nogil:

    cdef Tet_Patch tet_patch = tet_patches[item]

    cdef RayHitData hd = compute_tet_patch_hit(tet_patch.v, ray.org, ray.dir)

    # only count this is it's the closest hit
    if (hd.t < ray.tnear or hd.t > ray.Ng[0]):
        return

    if (hd.converged and 0 <= hd.u and 0 <= hd.v and hd.u + hd.v <= 1):

        # we have a hit, so update ray information
        ray.u = hd.u
        ray.v = hd.v
        ray.geomID = tet_patch.geomID
        ray.primID = item
        ray.Ng[0] = hd.t

        # sample the solution at the calculated point
        sample_tet10(tet_patches, ray)

    return
