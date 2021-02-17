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

    cdef int max_heapify(self, np.intp_t index) nogil except -1
    cdef int propagate_up(self, np.intp_t index) nogil except -1
    cdef int add(self, np.float64_t val) nogil except -1
    cdef int add_pid(self, np.float64_t val, np.int64_t pid) nogil except -1
    cdef int heap_append(self, np.float64_t val, np.int64_t ind) nogil except -1
    cdef np.float64_t extract_max(self) nogil except -1
    cdef int validate_heap(self) nogil except -1

cdef class NeighborList:
    cdef public np.float64_t[:] data
    cdef np.float64_t* data_ptr
    cdef public np.int64_t[:] pids
    cdef np.int64_t* pids_ptr
    cdef np.intp_t size
    cdef np.intp_t _max_size

    cdef int _update_memview(self) except -1
    cdef int _extend(self) nogil except -1
    cdef int add_pid(self, np.float64_t val, np.int64_t ind) nogil except -1
