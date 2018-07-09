"""
A cython implementation of the bounded priority queue

This is a priority queue that only keeps track of smallest k values that have
been added to it.

This priority queue is implemented with the configuration of having the largest
element at the beginning - this exploited to store nearest neighbour lists.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np

cimport cython

cdef class BoundedPriorityQueue:
    def __cinit__(self, np.intp_t max_elements, np.intp_t pids=0):
        self.max_elements = max_elements
        # mark invalid recently  values with -1
        self.heap = np.zeros(max_elements)-1
        self.heap_ptr = &(self.heap[0])
        # only allocate memory if we intend to store particle ID's
        if(pids == 1):
            self.pids = np.zeros(max_elements, dtype="int64")-1
            self.pids_ptr = &(self.pids[0])

        self.size = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int reset(self) nogil except -1:
        # utlity function useful to re-validate the pointers if we allocate
        # the heap to a memoryview
        self.pids_ptr = &(self.pids[0])
        self.heap_ptr = &(self.heap[0])
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int max_heapify(self, np.intp_t index) nogil except -1:
        # start at the top and move down making sure every value gets smaller
        # on the way, if not, we swap values
        cdef np.intp_t left = 2 * index + 1
        cdef np.intp_t right = 2 * index + 2
        cdef np.intp_t largest = index

        if left < self.size and self.heap_ptr[left] > self.heap_ptr[largest]:
            largest = left
        if right < self.size and self.heap_ptr[right] > self.heap_ptr[largest]:
            largest = right

        if largest != index:
            self.heap_ptr[index], self.heap_ptr[largest] = self.heap_ptr[largest], self.heap_ptr[index]
            self.max_heapify(largest)

        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int max_heapify_pid(self, np.intp_t index) nogil except -1:
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
            self.pids_ptr[index], self.pids_ptr[largest] = \
                self.pids_ptr[largest], self.pids_ptr[index]
            self.max_heapify_pid(largest)

        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int propagate_up(self, np.intp_t index) nogil except -1:
        # while a value is large than the parent, we swap with the parent and
        # move the value up
        while index != 0 and self.heap_ptr[(index - 1) // 2] < self.heap_ptr[index]:
            self.heap_ptr[index], self.heap_ptr[(index - 1) // 2] = \
                self.heap_ptr[(index - 1) // 2], self.heap_ptr[index]
            index = (index - 1) // 2

        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int propagate_up_pid(self, np.intp_t index) nogil except -1:
        while index != 0 and self.heap_ptr[(index - 1) // 2] < self.heap_ptr[index]:
            self.heap_ptr[index], self.heap_ptr[(index - 1) // 2] = self.heap_ptr[(index - 1) // 2], self.heap_ptr[index]
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
                self.heap_append(val)
        else:
            self.heap_append(val)
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int add_pid(self, np.float64_t val, np.int64_t ind) nogil except -1:
        cdef int i
        if self.size == self.max_elements:
            for i in range(self.size):
                if(ind == self.pids[i]):
                    return 0
            if val < self.heap_ptr[0]:
                self.extract_max_pid()
                self.heap_append_pid(val, ind)
        else:
            self.heap_append_pid(val, ind)
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int heap_append(self, np.float64_t val) nogil except -1:
        # we add an element to the end of the array, and tell the array the size
        # has increased. as this might not be the smallest element, we need to
        # propogate up
        self.heap_ptr[self.size] = val
        self.size += 1
        self.propagate_up(self.size - 1)
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int heap_append_pid(self, np.float64_t val, np.int64_t ind) nogil except -1:
        self.heap_ptr[self.size] = val
        self.pids_ptr[self.size] = ind
        self.size += 1
        self.propagate_up_pid(self.size - 1)
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef np.float64_t extract_max(self) nogil except -1:
        # we remove the maximum value from the array, and return what that value
        # was
        cdef np.float64_t maximum = self.heap_ptr[0]
        cdef np.float64_t val = self.heap_ptr[self.size-1]
        self.heap_ptr[self.size-1] = -1
        self.size -= 1
        if self.size > 0:
            self.heap_ptr[0] = val
            self.max_heapify(0)
        return maximum

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef np.float64_t extract_max_pid(self) nogil except -1:
        cdef np.float64_t maximum = self.heap_ptr[0]
        cdef np.float64_t val = self.heap_ptr[self.size-1]
        cdef np.int64_t ind = self.pids_ptr[self.size-1]
        self.heap_ptr[self.size-1] = -1
        self.pids_ptr[self.size-1] = -1
        self.size -= 1
        if self.size > 0:
            self.heap_ptr[0] = val
            self.pids_ptr[0] = ind
            self.max_heapify_pid(0)
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


def validate_pid():
    m = BoundedPriorityQueue(5, True)

    # Add elements to the queue
    elements = [0.1, 0.25, 1.33, 0.5, 3.2, 4.6, 2.0, 0.4, 4.0, .001]
    pids = [1,2,3,4,5,6,7,8,9,10]

    for el, pid in zip(elements, pids):
        m.add_pid(el, pid)

    m.extract_max_pid()
    m.extract_max_pid()
    m.extract_max_pid()

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
