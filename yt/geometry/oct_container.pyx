"""
Oct container

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

from libc.stdlib cimport malloc, free
cimport numpy as np
import numpy as np
from oct_container cimport Oct, OctAllocationContainer, OctreeContainer

cdef OctAllocationContainer *allocate_octs(
        int n_octs, OctAllocationContainer *prev):
    cdef OctAllocationContainer *n_cont
    cdef Oct *oct
    cdef int n, i, j, k
    n_cont = <OctAllocationContainer *> malloc(
        sizeof(OctAllocationContainer))
    n_cont.my_octs = <Oct *> malloc(sizeof(Oct) * n_octs)
    for n in range(n_octs):
        oct = &n_cont.my_octs[n]
        oct.parent = NULL
        oct.ind = oct.domain = -1
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    oct.children[i][j][k] = NULL
    if prev != NULL:
        prev.next = n_cont
    n_cont.next = NULL
    return n_cont

cdef void free_octs(
        OctAllocationContainer *first):
    cdef OctAllocationContainer *cur
    while first != NULL:
        cur = first
        free(first.my_octs)
        first = cur.next
        free(cur)

# Here is the strategy for RAMSES containers:
#   * Read each domain individually, creating *all* octs found in that domain
#     file, even if they reside on other CPUs.
#   * Only allocate octs that reside on >= domain
#   * For all octs, insert into tree, which may require traversing existing
#     octs
#   * Note that this doesn ot allow OctAllocationContainer to exactly be a
#     chunk, but it is close.  For IO chunking, we can theoretically examine
#     those octs that live inside a given allocator.

cdef class OctreeContainer:

    def __init__(self, domain_dimensions, domain_left_edge, domain_right_edge):
        self.nn[0], self.nn[1], self.nn[2] = domain_dimensions
        cdef int i, j, k, p
        p = 0
        self.nocts = 0 # Increment when initialized
        self.root_mesh = <Oct****> malloc(sizeof(void*) * self.nn[0])
        self.cont = allocate_octs(self.nn[0]*self.nn[1]*self.nn[2], NULL)
        for i in range(3):
            self.nn[i] = domain_dimensions[i]
        for i in range(self.nn[0]):
            self.root_mesh[i] = <Oct ***> malloc(sizeof(void*) * self.nn[1])
            for j in range(self.nn[1]):
                self.root_mesh[i][j] = <Oct **> malloc(sizeof(void*) * self.nn[2])
                for k in range(self.nn[2]):
                    self.root_mesh[i][j][k] = &self.cont.my_octs[p]
                    p += 1
        # We don't initialize the octs yet
        for i in range(3):
            self.DLE[i] = domain_left_edge[i]
            self.DRE[i] = domain_right_edge[i]

    def __dealloc__(self):
        free_octs(self.cont)
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                free(self.root_mesh[i][j])
            free(self.root_mesh[i])
        free(self.root_mesh)

    def add_ramses(self, int curdom, int curlevel, int ng,
                   np.ndarray[np.float64_t, ndim=2] pos,
                   np.ndarray[np.int64_t, ndim=1] index,
                   np.ndarray[np.int64_t, ndim=2] cpumap):
        cdef int level, no, p, i, j, k, oi, ind[3]
        cdef Oct* cur = self.root_mesh[0][0][0]
        cdef OctAllocationContainer *oa, *nextoa
        cdef np.float64_t pp[3], cp[3], dds[3]
        cdef int to_allocate = 0
        no = pos.shape[0]
        nextoa = self.cont
        while nextoa != NULL:
            oa = nextoa
            nextoa = oa.next
        to_allocate = ng
        if level > 0:
            oa = allocate_octs(to_allocate, oa)
        oi = 0
        # How do we bootstrap ourselves?
        for p in range(no):
            for i in range(3):
                pp[i] = pos[p, i]
                dds[i] = (self.DRE[i] + self.DLE[i])/self.nn[i]
                ind[i] = <int> (pp[i]/dds[i])
                cp[i] = (ind[i] + 0.5) * dds[i]
            cur = self.root_mesh[ind[0]][ind[1]][ind[2]]
            # Now we find the location we want
            # Note that RAMSES I think 1-indexes levels, but we don't.
            for level in range(curlevel):
                for i in range(3):
                    dds[i] = dds[i] / 2.0
                    if cp[i] > pp[i]: 
                        ind[i] = 0
                        cp[i] -= dds[i]/2.0
                    else:
                        ind[i] = 1
                        cp[i] += dds[i]/2.0
                # Check if it has not been allocated
                next = cur.children[ind[0]][ind[1]][ind[2]]
                if next == NULL:
                    cur.children[ind[0]][ind[1]][ind[2]] = &oa.my_octs[oi]
                    oi += 1
                    next = cur.children[ind[0]][ind[1]][ind[2]]
                    next.parent = cur
                cur = next
            cur.domain = curdom
            cur.ind = index[p]
            cur.local_ind = self.nocts
            self.nocts += 1
