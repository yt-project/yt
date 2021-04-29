"""
Bit array functions



"""


import numpy as np

cimport cython
cimport numpy as np


cdef inline void ba_set_value(np.uint8_t *buf, np.uint64_t ind,
                              np.uint8_t val) nogil:
    # This assumes 8 bit buffer
    if val > 0:
        buf[ind >> 3] |= (1 << (ind & 7))
    else:
        buf[ind >> 3] &= ~(1 << (ind & 7))

cdef inline np.uint8_t ba_get_value(np.uint8_t *buf, np.uint64_t ind) nogil:
    cdef np.uint8_t rv = (buf[ind >> 3] & (1 << (ind & 7)))
    if rv == 0: return 0
    return 1

cdef class bitarray:
    cdef np.uint8_t *buf
    cdef np.uint64_t size
    cdef np.uint64_t buf_size # Not exactly the same
    cdef public object ibuf

    cdef void _set_value(self, np.uint64_t ind, np.uint8_t val)
    cdef np.uint8_t _query_value(self, np.uint64_t ind)
    #cdef void set_range(self, np.uint64_t ind, np.uint64_t count, int val)
    #cdef int query_range(self, np.uint64_t ind, np.uint64_t count, int *val)
