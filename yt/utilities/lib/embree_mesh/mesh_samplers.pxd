cimport cython
cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr


cdef void sample_hex(void* userPtr,
                     rtcr.RTCRay& ray) nogil

cdef void sample_wedge(void* userPtr,
                       rtcr.RTCRay& ray) nogil

cdef void sample_tetra(void* userPtr,
                       rtcr.RTCRay& ray) nogil

cdef void sample_hex20(void* userPtr,
                       rtcr.RTCRay& ray) nogil

cdef void sample_tet10(void* userPtr,
                       rtcr.RTCRay& ray) nogil
