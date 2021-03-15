
# distutils: libraries = STD_LIBS
"""
An allocation container and memory pool



"""


cimport numpy as np

import numpy as np


cdef class ObjectPool:
    def __cinit__(self):
        """This class is *not* meant to be initialized directly, but instead
        through subclasses.  Those subclasses need to implement at a minimum
        the setting of itemsize, but can optionally also implement setup_objs
        and teardown_objs to either set default values or initialize additional
        pointers and values, and then free them."""
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
        """This returns an array view (if possible and implemented in a
        subclass) on the pool of objects specified by i."""
        return self._con_to_array(i)

    def __len__(self):
        return self.n_con

    def append(self, int n_objs, np.int64_t con_id = -1):
        """This allocates a new batch of n_objs objects, with the container id
        specified as con_id.  Return value is a view on them."""
        self.allocate_objs(n_objs, con_id)
        return self[self.n_con - 1]

    cdef void allocate_objs(self, int n_objs, np.int64_t con_id = -1) except *:
        cdef AllocationContainer *n_cont
        cdef AllocationContainer *prev
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

    cdef void setup_objs(self, void *obj, np.uint64_t count,
                         np.uint64_t offset, np.int64_t con_id):
        """This can be overridden in the subclass, where it can initialize any
        of the allocated objects."""
        pass

    cdef void teardown_objs(self, void *obj, np.uint64_t n, np.uint64_t offset,
                           np.int64_t con_id):
        # We assume that these are all allocated and have no sub-allocations
        """If overridden, additional behavior can be provided during teardown
        of allocations.  For instance, if you allocate some memory on each of
        the allocated objects."""
        if obj != NULL:
            free(obj)

    def to_arrays(self):
        rv = []
        cdef int i
        for i in range(self.n_con):
            rv.append(self._con_to_array(i))
        return rv

    def _con_to_array(self, int i):
        """This has to be implemented in a subclass, and should return an
        appropriate np.asarray() of a memoryview of self.my_objs."""
        raise NotImplementedError

cdef class BitmaskPool(ObjectPool):
    def __cinit__(self):
        """This is an example implementation of object pools for bitmasks
        (uint8) arrays.  It lets you reasonably quickly allocate a set of
        uint8 arrays which can be accessed and modified, and then virtually
        append to that."""
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
