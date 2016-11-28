"""
An allocation container and memory pool



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np
import numpy as np

cdef class ObjectPool:
    def __cinit__(self):
        self.containers = NULL
        self.n_con = 0
        self.itemsize = -1 # Never use the base class

    def __dealloc__(self):
        cdef int i
        cdef AllocationContainer *obj
        for i in range(self.n_con):
            obj = &self.containers[i]
            self.teardown_objs(obj.my_objs, obj.n, obj.offset, obj.con_id)
        if self.containers != NULL:
            free(self.containers)

    def __getitem__(self, int i):
        return self._con_to_array(i)

    def __len__(self):
        return self.n_con

    def append(self, int n_objs, np.int64_t con_id = -1):
        self.allocate_objs(n_objs, con_id)
        return self[self.n_con - 1]
        
    cdef void allocate_objs(self, int n_objs, np.int64_t con_id = -1) except *:
        cdef AllocationContainer *n_cont, *prev
        cdef int n, i, j, k
        cdef char *obj # char so we can do pointer math
        self.containers = <AllocationContainer*> realloc(
              self.containers, 
              sizeof(AllocationContainer) * (self.n_con + 1))
        n_cont = &self.containers[self.n_con]
        if self.n_con == 0:
            n_cont.offset = 0
        else:
            prev = &self.containers[self.n_con - 1]
            n_cont.offset = prev.offset + prev.n
        self.n_con += 1
        n_cont.my_objs = malloc(self.itemsize * n_objs)
        if n_cont.my_objs == NULL:
            raise MemoryError
        n_cont.n = n_objs
        n_cont.n_assigned = 0
        n_cont.con_id = con_id
        obj = <char*> n_cont.my_objs
        self.setup_objs(n_cont.my_objs, n_objs, n_cont.offset, n_cont.con_id)

    cdef void setup_objs(self, void *obj, np.uint64_t count, np.uint64_t
    offset, np.int64_t con_id):
        pass

    cdef void teardown_objs(self, void *obj, np.uint64_t n, np.uint64_t offset,
                           np.int64_t con_id):
        # We assume that these are all allocated and have no sub-allocations
        if obj != NULL:
            free(obj)

    def to_arrays(self):
        rv = []
        cdef int i
        for i in range(self.n_con):
            rv.append(self._con_to_array(i))
        return rv

    def _con_to_array(self, int i):
        raise NotImplementedError

cdef class BitmaskPool(ObjectPool):
    def __cinit__(self):
        # Base class will ALSO be called
        self.itemsize = sizeof(np.uint8_t)

    cdef void setup_objs(self, void *obj, np.uint64_t n, np.uint64_t offset,
                         np.int64_t con_id):
        cdef np.uint64_t i
        cdef np.uint8_t *mask = <np.uint8_t *> obj
        for i in range(n):
            mask[i] = 0

    def _con_to_array(self, int i):
        cdef AllocationContainer *obj = &self.containers[i]
        cdef np.uint8_t[:] o = <np.uint8_t[:obj.n]> (<np.uint8_t*> obj.my_objs)
        rv = np.asarray(o)
        return rv
