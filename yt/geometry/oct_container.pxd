"""
Oct definitions file

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

cimport numpy as np
from fp_utils cimport *
from selection_routines cimport SelectorObject
from oct_visitors cimport \
    OctVisitorData, oct_visitor_function, Oct
from libc.stdlib cimport bsearch, qsort

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
    cdef oct_visitor_function *fill_func
    cdef int partial_coverage
    cdef int level_offset
    cdef int nn[3]
    cdef np.uint8_t oref
    cdef np.float64_t DLE[3], DRE[3]
    cdef public np.int64_t nocts
    cdef public int num_domains
    cdef Oct *get(self, np.float64_t ppos[3], OctInfo *oinfo = ?)
    cdef int get_root(self, int ind[3], Oct **o)
    cdef Oct **neighbors(self, OctInfo *oinfo, np.int64_t *nneighbors)
    cdef void oct_bounds(self, Oct *, np.float64_t *, np.float64_t *)
    # This function must return the offset from global-to-local domains; i.e.,
    # OctAllocationContainer.offset if such a thing exists.
    cdef np.int64_t get_domain_offset(self, int domain_id)
    cdef void visit_all_octs(self, SelectorObject selector,
                        oct_visitor_function *func,
                        OctVisitorData *data)
    cdef Oct *next_root(self, int domain_id, int ind[3])
    cdef Oct *next_child(self, int domain_id, int ind[3], Oct *parent)
    cdef void setup_data(self, OctVisitorData *data, int domain_id = ?)

cdef class SparseOctreeContainer(OctreeContainer):
    cdef OctKey *root_nodes
    cdef void *tree_root
    cdef int num_root
    cdef int max_root
    cdef void key_to_ipos(self, np.int64_t key, np.int64_t pos[3])
    cdef np.int64_t ipos_to_key(self, int pos[3])

cdef class RAMSESOctreeContainer(SparseOctreeContainer):
    pass

cdef extern from "search.h" nogil:
    void *tsearch(const void *key, void **rootp,
                    int (*compar)(const void *, const void *))
    void *tfind(const void *key, const void **rootp,
                    int (*compar)(const void *, const void *))
    void *tdelete(const void *key, void **rootp,
                    int (*compar)(const void *, const void *))
