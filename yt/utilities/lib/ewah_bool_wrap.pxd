cimport numpy as np

cdef class BoolArrayCollection:
    cdef void* ewah_coll
    cdef void* ewah_keys

    cdef void _set(self, np.uint64_t i1, np.uint64_t i2)
    cdef bint _get(self, np.uint64_t i1, np.uint64_t i2)
