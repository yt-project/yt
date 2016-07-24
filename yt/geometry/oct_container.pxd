"""
Oct definitions file




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport cython
cimport numpy as np
from yt.utilities.lib.fp_utils cimport *
cimport oct_visitors
cimport selection_routines
from .oct_visitors cimport OctVisitor, Oct, cind
from libc.stdlib cimport bsearch, qsort, realloc, malloc, free
from libc.math cimport floor

cdef int ORDER_MAX

cdef struct OctKey:
    np.int64_t key
    Oct *node

cdef struct OctInfo:
    np.float64_t left_edge[3]
    np.float64_t dds[3]
    np.int64_t ipos[3]
    np.int32_t level

cdef struct OctAllocationContainer
cdef struct OctAllocationContainer:
    np.int64_t n
    np.int64_t n_assigned
    np.int64_t offset
    np.int64_t con_id
    OctAllocationContainer *next
    Oct *my_octs

cdef struct OctList

cdef struct OctList:
    OctList *next
    Oct *o

cdef OctList *OctList_append(OctList *list, Oct *o)
cdef int OctList_count(OctList *list)
cdef void OctList_delete(OctList *list)

cdef class OctreeContainer:
    cdef OctAllocationContainer *cont
    cdef OctAllocationContainer **domains
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
                  int max_level = ?)
    cdef int get_root(self, int ind[3], Oct **o)
    cdef Oct **neighbors(self, OctInfo *oinfo, np.int64_t *nneighbors,
                         Oct *o, bint periodicity[3])
    cdef void oct_bounds(self, Oct *, np.float64_t *, np.float64_t *)
    # This function must return the offset from global-to-local domains; i.e.,
    # OctAllocationContainer.offset if such a thing exists.
    cdef np.int64_t get_domain_offset(self, int domain_id)
    cdef void visit_all_octs(self,
                        selection_routines.SelectorObject selector,
                        OctVisitor visitor,
                        int vc = ?)
    cdef Oct *next_root(self, int domain_id, int ind[3])
    cdef Oct *next_child(self, int domain_id, int ind[3], Oct *parent)
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
    cdef np.int64_t ipos_to_key(self, int pos[3])

cdef class RAMSESOctreeContainer(SparseOctreeContainer):
    pass

cdef extern from "tsearch.h" nogil:
    void *tsearch(const void *key, void **rootp,
                    int (*compar)(const void *, const void *))
    void *tfind(const void *key, const void **rootp,
                    int (*compar)(const void *, const void *))
    void *tdelete(const void *key, void **rootp,
                    int (*compar)(const void *, const void *))
