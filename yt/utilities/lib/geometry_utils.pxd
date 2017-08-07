cimport numpy as np

cdef np.uint64_t position_to_morton(
    np.float64_t px, np.float64_t py, np.float64_t pz,
    np.float64_t dds[3], np.float64_t[:] dle, np.float64_t[:] dre) nogil

cdef int ORDER_MAX=20


