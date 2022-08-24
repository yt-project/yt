"""
A cython implementation of the bounded priority queue

This is a priority queue that only keeps track of smallest k values that have
been added to it.

This priority queue is implemented with the configuration of having the largest
element at the beginning - this exploited to store nearest neighbour lists.

"""


import numpy as np

cimport cython
cimport numpy as np
from cpython.mem cimport PyMem_Free, PyMem_Malloc, PyMem_Realloc


cdef class BoundedPriorityQueue:
    def __cinit__(self, np.intp_t max_elements, np.intp_t pids=0):
        self.max_elements = max_elements
        # mark invalid recently  values with -1
        self.heap = np.zeros(max_elements)-1
        self.heap_ptr = &(self.heap[0])
        # only allocate memory if we intend to store particle ID's
        self.use_pids = pids
        if pids == 1:
            self.pids = np.zeros(max_elements, dtype="int64")-1
            self.pids_ptr = &(self.pids[0])

        self.size = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int max_heapify(self, np.intp_t index) nogil except -1:
        cdef np.intp_t left = 2 * index + 1
        cdef np.intp_t right = 2 * index + 2
        cdef np.intp_t largest = index

        if left < self.size and self.heap_ptr[left] > self.heap_ptr[largest]:
            largest = left
        if right < self.size and self.heap_ptr[right] > self.heap_ptr[largest]:
            largest = right

        if largest != index:
            self.heap_ptr[index], self.heap_ptr[largest] = \
                self.heap_ptr[largest], self.heap_ptr[index]
            if self.use_pids:
                self.pids_ptr[index], self.pids_ptr[largest] = \
                    self.pids_ptr[largest], self.pids_ptr[index]

            self.max_heapify(largest)

        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int propagate_up(self, np.intp_t index) nogil except -1:
        while index != 0 and self.heap_ptr[(index - 1) // 2] < self.heap_ptr[index]:
            self.heap_ptr[index], self.heap_ptr[(index - 1) // 2] = self.heap_ptr[(index - 1) // 2], self.heap_ptr[index]
            if self.use_pids:
                self.pids_ptr[index], self.pids_ptr[(index - 1) // 2] = self.pids_ptr[(index - 1) // 2], self.pids_ptr[index]
            index = (index - 1) // 2

        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int add(self, np.float64_t val) nogil except -1:
        # if not at max size append, if at max size, only append if smaller than
        # the maximum value
        if self.size == self.max_elements:
            if val < self.heap_ptr[0]:
                self.extract_max()
                self.heap_append(val, -1)
        else:
            self.heap_append(val, -1)
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int add_pid(self, np.float64_t val, np.int64_t ind) nogil except -1:
        if self.size == self.max_elements:
            if val < self.heap_ptr[0]:
                self.extract_max()
                self.heap_append(val, ind)
        else:
            self.heap_append(val, ind)
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int heap_append(self, np.float64_t val, np.int64_t ind) nogil except -1:
        self.heap_ptr[self.size] = val
        if self.use_pids:
            self.pids_ptr[self.size] = ind
        self.size += 1
        self.propagate_up(self.size - 1)
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef np.float64_t extract_max(self) nogil except -1:
        cdef np.float64_t maximum = self.heap_ptr[0]
        cdef np.float64_t val
        cdef np.int64_t ind

        val = self.heap_ptr[self.size-1]
        self.heap_ptr[self.size-1] = -1

        if self.use_pids:
            ind = self.pids_ptr[self.size-1]
            self.pids_ptr[self.size-1] = -1

        self.size -= 1
        if self.size > 0:
            self.heap_ptr[0] = val
            if self.use_pids:
                self.pids_ptr[0] = ind
            self.max_heapify(0)
        return maximum

    cdef int validate_heap(self) nogil except -1:
        # this function loops through every element in the heap, if any children
        # are greater than their parents then we return zero, which is an error
        # as the heap condition is not satisfied
        cdef int i, index
        for i in range(self.size-1, -1, -1):
            index = i
            while index != 0:
                if self.heap_ptr[index] > self.heap_ptr[(index - 1) // 2]:
                    return 0
                index = (index - 1) // 2
        return 1

cdef class NeighborList:
    def __cinit__(self, np.intp_t init_size=32):
        self.size = 0
        self._max_size = init_size
        self.data_ptr = <np.float64_t*> PyMem_Malloc(
            self._max_size * sizeof(np.float64_t)
        )
        self.pids_ptr = <np.int64_t*> PyMem_Malloc(
            self._max_size * sizeof(np.int64_t)
        )
        self._update_memview()

    def __dealloc__(self):
        PyMem_Free(self.data_ptr)
        PyMem_Free(self.pids_ptr)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int _update_memview(self) except -1:
        self.data = <np.float64_t[:self._max_size]> self.data_ptr
        self.pids = <np.int64_t[:self._max_size]> self.pids_ptr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int _extend(self) nogil except -1:
        if self.size == self._max_size:
            self._max_size *= 2
            with gil:
                self.data_ptr = <np.float64_t*> PyMem_Realloc(
                    self.data_ptr,
                    self._max_size * sizeof(np.float64_t)
                )
                self.pids_ptr = <np.int64_t*> PyMem_Realloc(
                    self.pids_ptr,
                    self._max_size * sizeof(np.int64_t)
                )
                self._update_memview()
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int add_pid(self, np.float64_t val, np.int64_t ind) nogil except -1:
        self._extend()
        self.data_ptr[self.size] = val
        self.pids_ptr[self.size] = ind
        self.size += 1
        return 0

# these are test functions which are called from
# yt/utilities/lib/tests/test_nn.py
# they are stored here to allow easy interaction with functions not exposed at
# the python level
def validate_pid():
    m = BoundedPriorityQueue(5, True)

    # Add elements to the queue
    elements = [0.1, 0.25, 1.33, 0.5, 3.2, 4.6, 2.0, 0.4, 4.0, .001]
    pids = [1,2,3,4,5,6,7,8,9,10]

    for el, pid in zip(elements, pids):
        m.add_pid(el, pid)

    m.extract_max()
    m.extract_max()
    m.extract_max()

    return np.asarray(m.heap), np.asarray(m.pids)

def validate():
    m = BoundedPriorityQueue(5)

    # Add elements to the queue
    for el in [0.1, 0.25, 1.33, 0.5, 3.2, 4.6, 2.0, 0.4, 4.0, .001]:
        m.add(el)

    m.extract_max()
    m.extract_max()
    m.extract_max()

    return np.asarray(m.heap)

def validate_nblist():
    nblist = NeighborList(init_size=2)

    for i in range(4):
        nblist.add_pid(1.0, i)

    # Copy is necessary here. Without it, the allocated memory would be freed.
    # Leaving random data array.
    return np.asarray(nblist.data).copy(), np.asarray(nblist.pids).copy()
