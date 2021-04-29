cimport cython
cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_geometry as rtcg
cimport pyembree.rtcore_ray as rtcr

from .mesh_construction cimport Patch, Tet_Patch


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
