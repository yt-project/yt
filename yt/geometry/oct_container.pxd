"""
Oct definitions file




"""


cimport cython
cimport numpy as np
cimport oct_visitors
cimport selection_routines
from libc.math cimport floor
from libc.stdlib cimport bsearch, free, malloc, qsort, realloc

from yt.utilities.lib.allocation_container cimport AllocationContainer, ObjectPool
from yt.utilities.lib.fp_utils cimport *

from .oct_visitors cimport Oct, OctInfo, OctVisitor, cind


cdef int ORDER_MAX

cdef struct OctKey:
    np.int64_t key
    Oct *node
    # These next two are for particle sparse octrees.
    np.int64_t *indices
    np.int64_t pcount

cdef struct OctList

cdef struct OctList:
    OctList *next
    Oct *o

# NOTE: This object *has* to be the same size as the AllocationContainer
# object.  There's an assert in the __cinit__ function.
cdef struct OctAllocationContainer:
    np.uint64_t n
    np.uint64_t n_assigned
    np.uint64_t offset
    np.int64_t con_id # container id
    Oct *my_objs

cdef class OctObjectPool(ObjectPool):
    cdef inline OctAllocationContainer *get_cont(self, int i):
        return <OctAllocationContainer*> (&self.containers[i])

cdef OctList *OctList_append(OctList *list, Oct *o)
cdef int OctList_count(OctList *list)
cdef void OctList_delete(OctList *list)

cdef class OctreeContainer:
    cdef public OctObjectPool domains
    cdef Oct ****root_mesh
    cdef int partial_coverage
    cdef int level_offset
    cdef int nn[3]
    cdef np.uint8_t oref
    cdef np.float64_t DLE[3]
    cdef np.float64_t DRE[3]
    cdef public np.int64_t nocts
    cdef public int num_domains
    cdef Oct *get(self, np.float64_t ppos[3], OctInfo *oinfo = ?,
                  int max_level = ?) nogil
    cdef int get_root(self, int ind[3], Oct **o) nogil
    cdef Oct **neighbors(self, OctInfo *oinfo, np.int64_t *nneighbors,
                         Oct *o, bint periodicity[3])
    # This function must return the offset from global-to-local domains; i.e.,
    # AllocationContainer.offset if such a thing exists.
    cdef np.int64_t get_domain_offset(self, int domain_id)
    cdef void visit_all_octs(self,
                        selection_routines.SelectorObject selector,
                        OctVisitor visitor,
                        int vc = ?, np.int64_t *indices = ?)
    cdef Oct *next_root(self, int domain_id, int ind[3])
    cdef Oct *next_child(self, int domain_id, int ind[3], Oct *parent) except? NULL
    cdef void append_domain(self, np.int64_t domain_count)
    # The fill_style is the ordering, C or F, of the octs in the file.  "o"
    # corresponds to C, and "r" is for Fortran.
    cdef public object fill_style

cdef class SparseOctreeContainer(OctreeContainer):
    cdef OctKey *root_nodes
    cdef void *tree_root
    cdef int num_root
    cdef int max_root
    cdef void key_to_ipos(self, np.int64_t key, np.int64_t pos[3])
    cdef np.int64_t ipos_to_key(self, int pos[3]) nogil

cdef class RAMSESOctreeContainer(SparseOctreeContainer):
    pass

cdef extern from "tsearch.h" nogil:
    void *tsearch(const void *key, void **rootp,
                    int (*compar)(const void *, const void *))
    void *tfind(const void *key, const void **rootp,
                    int (*compar)(const void *, const void *))
    void *tdelete(const void *key, void **rootp,
                    int (*compar)(const void *, const void *))
