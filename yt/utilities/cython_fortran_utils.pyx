cimport numpy as np
import numpy as np
import cython
from libc.stdio cimport *
import struct

cdef INT32_SIZE = sizeof(np.int32_t)
cdef DOUBLE_SIZE = sizeof(np.float64_t)

cdef class FortranFile:
    def __cinit__(self, str fname):
        self.cfile = fopen(fname.encode('utf-8'), 'r')
        self._closed = False

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    cpdef void skip(self, INT64_t n=1):
        cdef INT32_t s1, s2, i

        if self._closed:
            raise ValueError("Read of closed file.")

        for i in range(n):
            fread(&s1, INT32_SIZE, 1, self.cfile)
            fseek(self.cfile, s1, SEEK_CUR)
            fread(&s2, INT32_SIZE, 1, self.cfile)

            if s1 != s2:
                raise IOError('Sizes do not agree in the header and footer for '
                              'this record - check header dtype')

    cdef INT64_t get_size(self, str dtype):
        if dtype == 'i':
            return 4
        elif dtype == 'd':
            return 8
        elif dtype == 'f':
            return 4
        elif dtype == 'l':
            return 8
        else:
            # Fallback to (slow) struct to compute the size
            return struct.calcsize(dtype)

    cpdef np.ndarray read_vector(self, str dtype):
        cdef INT32_t s1, s2, size
        cdef np.ndarray data

        if self._closed:
            raise ValueError("I/O operation on closed file.")

        size = self.get_size(dtype)

        fread(&s1, INT32_SIZE, 1, self.cfile)

        # Check record is compatible with data type
        if s1 % size != 0:
            raise ValueError('Size obtained (%s) does not match with the expected '
                             'size (%s) of multi-item record', s1, size)

        data = np.empty(s1 // size, dtype=dtype)
        fread(<void *>data.data, size, s1 // size, self.cfile)
        fread(&s2, INT32_SIZE, 1, self.cfile)

        if s1 != s2:
            raise IOError('Sizes do not agree in the header and footer for '
                          'this record - check header dtype')

        return data

    cpdef INT32_t read_int(self):
        cdef INT32_t s1, s2
        cdef INT32_t data

        if self._closed:
            raise ValueError("I/O operation on closed file.")

        fread(&s1, INT32_SIZE, 1, self.cfile)

        if s1 % INT32_SIZE != 0:
            raise ValueError('Size obtained (%s) does not match with the expected '
                             'size (%s) of multi-item record', s1, INT32_SIZE)

        fread(&data, INT32_SIZE, s1 // INT32_SIZE, self.cfile)
        fread(&s2, INT32_SIZE, 1, self.cfile)

        if s1 != s2:
            raise IOError('Sizes do not agree in the header and footer for '
                          'this record - check header dtype')

        return data

    cpdef dict read_attrs(self, object attrs):
        cdef str dtype
        cdef int n
        cdef dict data
        cdef key
        cdef np.ndarray tmp

        if self._closed:
            raise ValueError("I/O operation on closed file.")

        data = {}

        for key, n, dtype in attrs:
            if n == 1:
                data[key] = self.read_vector(dtype)[0]
            else:
                tmp = self.read_vector(dtype)
                if type(key) == tuple:
                    # There are multiple keys
                    for ikey in range(n):
                        data[key[ikey]] = tmp[ikey]
                else:
                    data[key] = tmp

        return data

    cpdef INT64_t tell(self):
        cdef INT64_t pos
        if self._closed:
            raise ValueError("I/O operation on closed file.")

        pos = ftell(self.cfile)
        return pos

    cpdef void seek(self, INT64_t pos, INT64_t whence=SEEK_SET):
        if self._closed:
            raise ValueError("I/O operation on closed file.")
        fseek(self.cfile, pos, whence)

    cpdef void close(self):
        if self._closed:
            return
        fclose(self.cfile)
        self._closed = True

    def __dealloc__(self):
        if not self._closed:
            self.close()
