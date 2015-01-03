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
from libc.stdlib cimport malloc, free

cdef class bitarray:

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __init__(self, size = -1, arr = None):
        cdef np.uint64_t i
        if size == -1 and arr is None:
            raise RuntimeError
        elif size == -1:
            size = arr.size
        elif size != -1 and arr is not None:
            if size != arr.size:
                raise RuntimeError
        self.buf_size = (size >> 3)
        if (size & 7) != 0:
            # We need an extra one if we've got any lingering bits
            self.buf_size += 1
        cdef np.ndarray[np.uint8_t] ibuf_t
        ibuf_t = self.ibuf = np.zeros(self.buf_size, "uint8")
        self.buf = <np.uint8_t *> ibuf_t.data
        self.size = size
        if arr is not None:
            self.set_from_array(arr)
        else:
            for i in range(self.buf_size):
                self.buf[i] = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def set_from_array(self, np.ndarray[np.uint8_t, cast=True] arr):
        cdef np.uint64_t i, j, elem
        cdef np.uint8_t *btemp = self.buf
        arr = np.ascontiguousarray(arr)
        j = 0
        for i in range(self.size):
            btemp[i >> 3] = btemp[i >> 3] | (arr[i] << (7-j))
            j += 1
            if j == 8:
                j = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def as_bool_array(self):
        cdef np.uint64_t i, j, elem
        cdef np.uint8_t *btemp = self.buf
        cdef np.ndarray[np.uint8_t, ndim=1] output
        output = np.zeros(self.size, "uint8")
        j = 0
        for i in range(self.size):
            output[i] = (btemp[i >> 3] >> (7 - j)) & 1
            j += 1
            if j == 8:
                j = 0
        return output.astype("bool")

    cdef void _set_value(self, np.uint64_t ind, np.uint8_t val):
        ba_set_value(self.buf, ind, val)

    def set_value(self, np.uint64_t ind, np.uint8_t val):
        ba_set_value(self.buf, ind, val)

    cdef np.uint8_t _query_value(self, np.uint64_t ind):
        return ba_get_value(self.buf, ind)

    def query_value(self, np.uint64_t ind):
        return ba_get_value(self.buf, ind)
