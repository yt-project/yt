cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr
cimport cython

cdef void sample_hex(void* userPtr,
                     rtcr.RTCRay& ray) nogil

cdef void sample_tetra(void* userPtr,
                       rtcr.RTCRay& ray) nogil

cdef double sample_hex_at_real_point(double* vertices,
                                     double* field_values,
                                     double* physical_x) nogil
