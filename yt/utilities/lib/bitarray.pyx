# distutils: libraries = STD_LIBS
"""
Bit array functions



"""


import numpy as np

cimport cython
cimport numpy as np


cdef class bitarray:

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __cinit__(self, np.int64_t size = -1,
                  np.ndarray[np.uint8_t, ndim=1, cast=True] arr = None):
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
        >>> print(a & b)
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
        cdef np.uint8_t bitmask = 255
        if (size & 7) != 0:
            # We need an extra one if we've got any lingering bits
            self.buf_size += 1
            bitmask = 0
            for i in range(size & 7):
                bitmask |= (1<<i)
        self.final_bitmask = bitmask
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
    def set_from_array(self, np.ndarray[np.uint8_t, cast=True] arr not None):
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
        >>> print(a.set_value(2, 1))

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
        >>> print(a.query_value(2))

        """
        return ba_get_value(self.buf, ind)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void _set_range(self, np.uint64_t start, np.uint64_t stop, np.uint8_t val):
        ba_set_range(self.buf, start, stop, val)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def set_range(self, np.uint64_t start, np.uint64_t stop, np.uint8_t val):
        r"""Set a range of values to on/off.  Uses slice-style indexing.

        No return value.

        Parameters
        ----------
        start : int
            The starting component of a slice.
        stop : int
            The ending component of a slice.
        val : bool or uint8_t
            What to set the range to

        Examples
        --------

        >>> arr_in = np.array([True, True, False, True, True, False])
        >>> a = ba.bitarray(arr = arr_in)
        >>> a.set_range(0, 3, 0)

        """
        ba_set_range(self.buf, start, stop, val)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.uint64_t _count(self):
        cdef np.uint64_t count = 0
        cdef np.uint64_t i
        self.buf[self.buf_size - 1] &= self.final_bitmask
        for i in range(self.buf_size):
            count += _num_set_bits(self.buf[i])
        return count


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count(self):
        r"""Count the number of values set in the array.

        Parameters
        ----------

        Examples
        --------

        >>> arr_in = np.array([True, True, False, True, True, False])
        >>> a = ba.bitarray(arr = arr_in)
        >>> a.count()

        """
        return self._count()

    cdef bitarray _logical_and(self, bitarray other, bitarray result = None):
        # Create a place to put it.  Note that we might have trailing values,
        # we actually need to reset the ending set.
        if other.size != self.size:
            raise IndexError
        if result is None:
            result = bitarray(self.size)
        for i in range(self.buf_size):
            result.buf[i] = other.buf[i] & self.buf[i]
        result.buf[self.buf_size - 1] &= self.final_bitmask
        return result

    def logical_and(self, bitarray other, bitarray result = None):
        return self._logical_and(other, result)

    def __and__(self, bitarray other):
        # Wrap it directly here.
        return self.logical_and(other)

    def __iand__(self, bitarray other):
        rv = self.logical_and(other, self)
        return rv

    cdef bitarray _logical_or(self, bitarray other, bitarray result = None):
        if other.size != self.size:
            raise IndexError
        if result is None:
            result = bitarray(self.size)
        for i in range(self.buf_size):
            result.buf[i] = other.buf[i] | self.buf[i]
        result.buf[self.buf_size - 1] &= self.final_bitmask
        return result

    def logical_or(self, bitarray other, bitarray result = None):
        return self._logical_or(other, result)

    def __or__(self, bitarray other):
        return self.logical_or(other)

    def __ior__(self, bitarray other):
        return self.logical_or(other, self)

    cdef bitarray _logical_xor(self, bitarray other, bitarray result = None):
        if other.size != self.size:
            raise IndexError
        if result is None:
            result = bitarray(self.size)
        for i in range(self.buf_size):
            result.buf[i] = other.buf[i] ^ self.buf[i]
        result.buf[self.buf_size - 1] &= self.final_bitmask
        return result

    def logical_xor(self, bitarray other, bitarray result = None):
        return self._logical_xor(other, result)

    def __xor__(self, bitarray other):
        return self.logical_xor(other)

    def __ixor__(self, bitarray other):
        return self.logical_xor(other, self)
