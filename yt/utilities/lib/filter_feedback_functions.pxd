cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr
from pyembree.rtcore cimport Vec3f
cimport cython

cdef double sample_surface_hex(void* userPtr,
                               rtcr.RTCRay& ray)

cdef double get_value_trilinear(void* userPtr,
                                rtcr.RTCRay& ray)

cdef void maximum_intensity(void* userPtr, 
                            rtcr.RTCRay& ray)

cdef void sample_surface(void* userPtr, rtcr.RTCRay& ray)
