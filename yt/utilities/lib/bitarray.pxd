"""
Bit array functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython

cdef inline void ba_set_value(np.uint8_t *buf, np.uint64_t ind,
                              np.uint8_t val) nogil:
    if val > 0: val = 1
    buf[ind >> 3] |= (val << (ind & 7))

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
