cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr
cimport pyembree.rtcore_geometry as rtcg
from yt.utilities.lib.mesh_construction cimport Patch
cimport cython

cdef void patchIntersectFunc(Patch* patches,
                             rtcr.RTCRay& ray,
                             size_t item) nogil

cdef void patchBoundsFunc(Patch* patches,
                          size_t item,
                          rtcg.RTCBounds* bounds_o) nogil
