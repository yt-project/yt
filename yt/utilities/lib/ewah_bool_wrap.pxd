cimport numpy as np

cdef class BoolArrayCollection:
    cdef void* ewah_coll
    cdef void* ewah_keys
    cdef void* ewah_refn
    cdef void* ewah_coar

    cdef bint _richcmp(self, BoolArrayCollection solf, int op)
    cdef void _set(self, np.uint64_t i1, np.uint64_t i2=*)
    cdef void _set_coarse(self, np.uint64_t i1)
    cdef void _set_refined(self, np.uint64_t i1, np.uint64_t i2)
    cdef void _set_map(self, np.uint64_t i1, np.uint64_t i2)
    cdef void _set_refn(self, np.uint64_t i1)
    cdef bint _get(self, np.uint64_t i1, np.uint64_t i2=*)
    cdef bint _get_coarse(self, np.uint64_t i1)
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
    cdef int total
    cdef void* entries
    cdef void _set(self, np.uint64_t ind)
    cdef void _fill(self, np.uint8_t[:] mask)
    cdef void _fill_ewah(self, BoolArrayCollection mm)
    cdef void _reset(self)
    cdef to_array(self)
    cdef void _remove_duplicates(self)
    cdef void _prune(self)

cdef class SparseUnorderedRefinedBitmask:
    cdef int total
    cdef void* entries
    # cdef void* entries1
    # cdef void* entries2
    cdef void _set(self, np.uint64_t ind1, np.uint64_t ind2)
    cdef void _fill(self, np.uint8_t[:] mask1, np.uint8_t[:])
    cdef void _fill_ewah(self, BoolArrayCollection mm)
    cdef void _reset(self)
    cdef to_array(self)
    cdef void _remove_duplicates(self)
    cdef void _prune(self)
