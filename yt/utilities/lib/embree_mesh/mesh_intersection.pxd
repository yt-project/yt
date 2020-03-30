cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr
cimport pyembree.rtcore_geometry as rtcg
from yt.utilities.lib.mesh_construction cimport Patch, Tet_Patch
cimport cython

cdef void patchIntersectFunc(Patch* patches,
                             rtcr.RTCRay& ray,
                             size_t item) nogil

cdef void patchBoundsFunc(Patch* patches,
                          size_t item,
                          rtcg.RTCBounds* bounds_o) nogil

cdef void tet_patchIntersectFunc(Tet_Patch* tet_patches,
                             rtcr.RTCRay& ray,
                             size_t item) nogil

cdef void tet_patchBoundsFunc(Tet_Patch* tet_patches,
                          size_t item,
                          rtcg.RTCBounds* bounds_o) nogil
