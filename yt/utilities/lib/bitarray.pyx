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
        r"""This is a bitarray, which flips individual bits to on/off inside a
        uint8 container array.

        By encoding on/off inside each bit in a uint8 array, we can compress
        boolean information down by up to a factor of 8.  Either an input array
        or a size must be provided.

        Parameters
        ----------
        size : int
            The size we should pre-allocate.
        arr : array-like
            An input array to turn into a bitarray.

        Examples
        --------

        >>> arr_in1 = np.array([True, True, False])
        >>> arr_in2 = np.array([False, True, True])
        >>> a = ba.bitarray(arr = arr_in1)
        >>> b = ba.bitarray(arr = arr_in2)
        >>> print a & b
        >>> print (a & b).as_bool_array()

        """
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
        r"""Given an array that is either uint8_t or boolean, set the values of
        this array to match it.

        Parameters
        ----------
        arr : array, castable to uint8
            The array we set from.
        """
        cdef np.uint64_t i, j
        cdef np.uint8_t *btemp = self.buf
        arr = np.ascontiguousarray(arr)
        j = 0
        for i in range(self.size):
            btemp[i >> 3] = btemp[i >> 3] | (arr[i] << (j))
            j += 1
            if j == 8:
                j = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def as_bool_array(self):
        r"""Return a copy of this array, as a boolean array.

        All of the values encoded in this bitarray are expanded into boolean
        values in a new array and returned.

        Returns
        -------
        arr : numpy array of type bool
            The uint8 values expanded into boolean values

        """
        cdef np.uint64_t i, j
        cdef np.uint8_t *btemp = self.buf
        cdef np.ndarray[np.uint8_t, ndim=1] output
        output = np.zeros(self.size, "uint8")
        j = 0
        for i in range(self.size):
            output[i] = (btemp[i >> 3] >> (j)) & 1
            j += 1
            if j == 8:
                j = 0
        return output.astype("bool")

    cdef void _set_value(self, np.uint64_t ind, np.uint8_t val):
        ba_set_value(self.buf, ind, val)

    def set_value(self, np.uint64_t ind, np.uint8_t val):
        r"""Set the on/off value of a given bit.

        Modify the value encoded in a given index.

        Parameters
        ----------
        ind : int
            The index to query in the bitarray.
        val : bool or uint8_t
            What to set the index to

        Examples
        --------

        >>> arr_in = np.array([True, True, False])
        >>> a = ba.bitarray(arr = arr_in)
        >>> print a.set_value(2, 1)

        """
        ba_set_value(self.buf, ind, val)

    cdef np.uint8_t _query_value(self, np.uint64_t ind):
        return ba_get_value(self.buf, ind)

    def query_value(self, np.uint64_t ind):
        r"""Query the on/off value of a given bit.

        Return the value encoded in a given index.

        Parameters
        ----------
        ind : int
            The index to query in the bitarray.

        Examples
        --------

        >>> arr_in = np.array([True, True, False])
        >>> a = ba.bitarray(arr = arr_in)
        >>> print a.query_value(2)

        """
        return ba_get_value(self.buf, ind)
