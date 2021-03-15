"""
An allocation container and memory pool



"""


cimport numpy as np
from libc.stdlib cimport free, malloc, realloc


cdef struct AllocationContainer:
    np.uint64_t n
    np.uint64_t n_assigned
    np.uint64_t offset
    np.int64_t con_id # container id
    void *my_objs

cdef class ObjectPool:
    cdef public np.uint64_t itemsize
    cdef np.uint64_t n_con
    cdef AllocationContainer* containers
    cdef void allocate_objs(self, int n_objs, np.int64_t con_id = ?) except *
    cdef void setup_objs(self, void *obj, np.uint64_t count,
                         np.uint64_t offset, np.int64_t con_id)
    cdef void teardown_objs(self, void *obj, np.uint64_t n, np.uint64_t offset,
                           np.int64_t con_id)
