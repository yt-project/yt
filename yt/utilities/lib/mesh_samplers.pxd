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

cdef double sample_tetra_at_real_point(double* vertices,
                                       double* field_values,
                                       double* physical_x) nogil

cdef void hex_real_to_mapped(double* mapped_x,
                             double* vertices,
                             double* physical_x) nogil

cdef double sample_hex_at_unit_point(double* coord, double* vals) nogil

cdef int hex_check_inside(double* mapped_coord) nogil

cdef void tetra_real_to_mapped(double* mapped_coord,
                               double* vertices,
                               double* physical_coord) nogil

cdef double sample_tetra_at_unit_point(double* coord, double* vals) nogil

cdef int tetra_check_inside(double* mapped_coord) nogil
