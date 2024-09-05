"""
A cython implementation of the bounded priority queue

This is a priority queue that only keeps track of smallest k values that have
been added to it.


"""

import numpy as np

cimport numpy as np


cdef class BoundedPriorityQueue:
    cdef public np.float64_t[:] heap
    cdef np.float64_t* heap_ptr
    cdef public np.int64_t[:] pids
    cdef np.int64_t* pids_ptr
    cdef int use_pids

    cdef np.intp_t size
    cdef np.intp_t max_elements

    cdef int max_heapify(self, np.intp_t index) except -1 nogil
    cdef int propagate_up(self, np.intp_t index) except -1 nogil
    cdef int add(self, np.float64_t val) except -1 nogil
    cdef int add_pid(self, np.float64_t val, np.int64_t pid) except -1 nogil
    cdef int heap_append(self, np.float64_t val, np.int64_t ind) except -1 nogil
    cdef np.float64_t extract_max(self) except -1 nogil
    cdef int validate_heap(self) except -1 nogil

cdef class NeighborList:
    cdef public np.float64_t[:] data
    cdef np.float64_t* data_ptr
    cdef public np.int64_t[:] pids
    cdef np.int64_t* pids_ptr
    cdef np.intp_t size
    cdef np.intp_t _max_size

    cdef int _update_memview(self) except -1
    cdef int _extend(self) except -1 nogil
    cdef int add_pid(self, np.float64_t val, np.int64_t ind) except -1 nogil
