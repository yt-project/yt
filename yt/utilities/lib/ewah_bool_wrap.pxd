cimport numpy as np

cdef class BoolArrayCollection:
    cdef void* ewah_coll
    cdef void* ewah_keys
    cdef void* ewah_refn
    cdef void* ewah_coar

    cdef void _set(self, np.uint64_t i1, np.uint64_t i2=*)
    cdef void _set_map(self, np.uint64_t i1, np.uint64_t i2)
    cdef void _set_refn(self, np.uint64_t i1)
    cdef bint _get(self, np.uint64_t i1, np.uint64_t i2=*)
    cdef bint _contains(self, np.uint64_t i)
    cdef bint _isref(self, np.uint64_t i)
    cdef void _ewah_coarse(self)
    cdef int _count_total(self)
    cdef int _count_refined(self)
    cdef int _count_coarse(self)
    cdef void _append(self, BoolArrayCollection solf)
    cdef bytes _dumps(self)
    cdef void _loads(self, bytes s)

cdef class SparseUnorderedBitmask:
    cdef void* entries
    cdef void _set(self, np.uint64_t ind)
    cdef void _fill(self, np.uint8_t[:] mask)
