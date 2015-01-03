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

cdef class bitarray:
    cdef np.uint8_t *buf
    cdef np.uint64_t size
    cdef np.uint64_t buf_size # Not exactly the same
    cdef public object ibuf

    #cpdef int set_value(self, np.uint64_t ind, int val)
    #cpdef int query_value(self, np.uint64_t ind)
    #cpdef int set_range(self, np.uint64_t ind, np.uint64_t count, int val)
    #cpdef int query_range(self, np.uint64_t ind, np.uint64_t count, int *val)
