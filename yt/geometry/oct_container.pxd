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
from selection_routines cimport SelectorObject, \
    OctVisitorData, oct_visitor_function
from oct_visitors cimport *
from libc.stdlib cimport bsearch, qsort

cdef int ORDER_MAX

cdef struct Oct
cdef struct Oct:
    np.int64_t file_ind     # index with respect to the order in which it was
                            # added
    np.int64_t domain_ind   # index within the global set of domains
                            # note that moving to a local index will require
                            # moving to split-up masks, which is part of a
                            # bigger refactor
    np.int64_t domain       # (opt) addl int index
    np.int64_t pos[3]       # position in ints
    np.int8_t level
    Oct *children[2][2][2]
    Oct *parent

cdef struct OctKey:
    np.int64_t key
    Oct *node

cdef struct OctInfo:
    np.float64_t left_edge[3]
    np.float64_t dds[3]

cdef struct OctAllocationContainer
cdef struct OctAllocationContainer:
    np.int64_t n
    np.int64_t n_assigned
    np.int64_t offset
    OctAllocationContainer *next
    Oct *my_octs

cdef class OctreeContainer:
    cdef OctAllocationContainer *cont
    cdef Oct ****root_mesh
    cdef int nn[3]
    cdef np.float64_t DLE[3], DRE[3]
    cdef public int nocts
    cdef public int max_domain
    cdef Oct *get(self, np.float64_t ppos[3], OctInfo *oinfo = ?)
    cdef int get_root(self, int ind[3], Oct **o)
    cdef void neighbors(self, Oct *, Oct **)
    cdef void oct_bounds(self, Oct *, np.float64_t *, np.float64_t *)
    # This function must return the offset from global-to-local domains; i.e.,
    # OctAllocationContainer.offset if such a thing exists.
    cdef np.int64_t get_domain_offset(self, int domain_id)
    cdef void visit_all_octs(self, SelectorObject selector,
                        oct_visitor_function *func,
                        OctVisitorData *data)

cdef class RAMSESOctreeContainer(OctreeContainer):
    cdef OctAllocationContainer **domains
    cdef OctKey *root_nodes
    cdef void *tree_root
    cdef int num_root
    cdef int max_root
    cdef Oct *next_root(self, int domain_id, int ind[3])
    cdef Oct *next_child(self, int domain_id, int ind[3], Oct *parent)

cdef extern from "search.h" nogil:
    void *tsearch(const void *key, void **rootp,
                    int (*compar)(const void *, const void *))
    void *tfind(const void *key, const void **rootp,
                    int (*compar)(const void *, const void *))
    void *tdelete(const void *key, void **rootp,
                    int (*compar)(const void *, const void *))

cdef inline np.int64_t oct_key(Oct *o):
    cdef int i
    if o.level != 0: return -1
    cdef np.int64_t key = 0
    for i in range(3):
        key |= (o.pos[i] << 20 * (2 - i))
    return key
