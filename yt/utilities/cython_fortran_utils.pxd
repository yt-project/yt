from libc.stdio cimport FILE
cimport numpy as np

cdef class FortranFile:
    cdef FILE* cfile

    cpdef void skip(self, int n=*)
    cdef int get_size(self, str dtype)
    cpdef np.ndarray read_vector(self, str dtype)
    cpdef int read_int(self)
    cpdef dict read_attrs(self, object attrs)
    cpdef long tell(self)
    cpdef void seek(self, int pos, int whence=*)
    cpdef void close(self)
