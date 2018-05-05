cimport numpy as np
import numpy as np
import cython
from libc.stdio cimport *
import struct

cdef INT_SIZE = sizeof(int)
cdef DOUBLE_SIZE = sizeof(double)


cdef class FortranFile:
    def __init__(self, str fname):
        self.cfile = fopen(fname.encode('utf-8'), 'r')

    cpdef void skip(self, int n=1):
        cdef int s1, s2, i

        for i in range(n):
            fread(&s1, INT_SIZE, 1, self.cfile)
            fseek(self.cfile, s1, SEEK_CUR)
            fread(&s2, INT_SIZE, 1, self.cfile)

    cdef int get_size(self, str dtype):
        if dtype == 'i':
            return 4
        elif dtype == 'd':
            return 8
        elif dtype == 'f':
            return 4
        else:
            # Fallback to struct to compute the size
            return struct.calcsize(dtype)

    cpdef np.ndarray read_vector(self, str dtype):
        cdef int s1, s2, size
        cdef np.ndarray data

        size = self.get_size(dtype)

        fread(&s1, INT_SIZE, 1, self.cfile)
        data = np.empty(s1 // size, dtype=dtype)
        fread(<void *>data.data, size, s1 // size, self.cfile)
        fread(&s2, INT_SIZE, 1, self.cfile)

        return data

    cpdef int read_int(self):
        cdef int s1, s2, size=sizeof(int)
        cdef int data

        fread(&s1, INT_SIZE, 1, self.cfile)
        fread(&data, size, s1 // size, self.cfile)
        fread(&s2, INT_SIZE, 1, self.cfile)

        return data

    cpdef dict read_attrs(self, object attrs):
        cdef str dtype
        cdef int n
        cdef dict data
        cdef key
        cdef np.ndarray tmp

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

    cpdef long tell(self):
        cdef long pos
        pos = ftell(self.cfile)
        return pos

    cpdef void seek(self, int pos, int whence=SEEK_SET):
        fseek(self.cfile, pos, whence)

    cpdef void close(self):
        fclose(self.cfile)
