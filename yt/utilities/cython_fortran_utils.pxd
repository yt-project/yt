from libc.stdio cimport FILE
cimport numpy as np

ctypedef np.int32_t INT32_t
ctypedef np.int64_t INT64_t
ctypedef np.float64_t DOUBLE_t

cdef class FortranFile:
    cdef FILE* cfile
    cdef bint _opened

    cpdef void skip(self, INT64_t n=*)
    cdef INT64_t get_size(self, str dtype)
    cpdef np.ndarray read_vector(self, str dtype)
    cpdef INT32_t read_int(self)
    cpdef dict read_attrs(self, object attrs)
    cpdef INT64_t tell(self)
    cpdef void seek(self, INT64_t pos, INT64_t whence=*)
    cpdef void close(self)
