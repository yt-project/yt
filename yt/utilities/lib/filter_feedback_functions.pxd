cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr
from pyembree.rtcore cimport Vec3f
cimport cython

cdef void sample_hex(void* userPtr,
                     rtcr.RTCRay& ray) nogil
