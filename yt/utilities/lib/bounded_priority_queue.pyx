"""
A cython implementation of the bounded priority queue

This is a priority queue that only keeps track of smallest k values that have
been added to it.


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
    def __cinit__(self, np.intp_t max_elements):
        self.max_elements = max_elements
        # mark invalid recently  values with -1
        self.heap = np.zeros(max_elements)-1
        self.heap_ptr = &(self.heap[0])
        self.size = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int max_heapify(self, np.intp_t index) nogil except -1:
        # this uses the traditional mapping of a heap onto an array with left
        # and rights
        # essentially just start at the top and move down making sure everything
        # gets smaller on the way
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
    cdef int propagate_up(self, np.intp_t index) nogil except -1:
        # while we are maller than a parent, we move the value up through the
        # array by swapping elements
        while index != 0 and self.heap_ptr[(index - 1) // 2] < self.heap_ptr[index]:
            self.heap_ptr[index], self.heap_ptr[(index - 1) // 2] = self.heap_ptr[(index - 1) // 2], self.heap_ptr[index]
            index = (index - 1) // 2

        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef int add(self, np.float64_t val) nogil except -1:
        # if not at max size append, if at max size, only append if smaller than
        # the maximum valuei
        # we append to the top and then we need to sift the value down to
        # validate the heap
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
    cdef np.float64_t extract_max(self) nogil except -1:
        cdef np.float64_t maximum = self.heap_ptr[0]
        cdef np.float64_t val = self.heap_ptr[self.size-1]
        self.heap_ptr[self.size-1] = -1
        self.size -= 1
        if self.size > 0:
            self.heap_ptr[0] = val
            self.max_heapify(0)
        return maximum

    # this function loops through every element in theheap, if any children are
    # greater than their parents then we return zero, which is an error and the
    # heap condition is not satisfied
    cdef int validate_heap(self) nogil except -1:
        cdef int i, index
        for i in range(self.size-1, -1, -1):
            index = i
            while index != 0:
                if self.heap_ptr[index] > self.heap_ptr[(index - 1) // 2]:
                    return 0
                index = (index - 1) // 2
        return 1

def validate():
    m = BoundedPriorityQueue(5)

    print("Initial heap:", np.asarray(m.heap))

    # Add elements to the queue
    for el in [0.1, 0.25, 1.33, 0.5, 3.2, 4.6, 2.0, 0.4, 4.0, .001]:
        m.add(el)
        print(np.asarray(m.heap))

    print(np.asarray(m.heap))

    print("Extract maximum:", m.extract_max())
    print("Extract maximum:", m.extract_max())
    print("Extract maximum:", m.extract_max())
    print("Extract maximum:", m.extract_max())

    print(np.asarray(m.heap))
