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
from selection_routines cimport SelectorObject
cimport cython

cdef extern from "stdlib.h":
    # NOTE that size_t might not be int
    void *alloca(int)

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
        oct.local_ind = -1
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

    def __init__(self, domain_dimensions, domain_left_edge, domain_right_edge,
                 int initial_allocation = 0):
        self.nn[0], self.nn[1], self.nn[2] = domain_dimensions
        cdef int i, j, k, p
        self.max_domain = -1
        p = 0
        self.nocts = 0 # Increment when initialized
        self.root_mesh = <Oct****> malloc(sizeof(void*) * self.nn[0])
        if initial_allocation == 0:
            self.cont = allocate_octs(self.nn[0]*self.nn[1]*self.nn[2], NULL)
        else:
            self.cont = allocate_octs(initial_allocation, NULL)
            for i in range(initial_allocation):
                self.cont.my_octs[i].local_ind = i
        for i in range(3):
            self.nn[i] = domain_dimensions[i]
        for i in range(self.nn[0]):
            self.root_mesh[i] = <Oct ***> malloc(sizeof(void*) * self.nn[1])
            for j in range(self.nn[1]):
                self.root_mesh[i][j] = <Oct **> malloc(sizeof(void*) * self.nn[2])
                for k in range(self.nn[2]):
                    self.root_mesh[i][j][k] = &self.cont.my_octs[self.nocts]
                    self.root_mesh[i][j][k].pos[0] = i
                    self.root_mesh[i][j][k].pos[1] = j
                    self.root_mesh[i][j][k].pos[2] = k
                    self.nocts += 1
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

    def __iter__(self):
        cdef OctAllocationContainer *cur = self.cont
        cdef Oct *this
        cdef int i
        while cur != NULL:
            for i in range(cur.n):
                this = &cur.my_octs[i]
                yield (this.ind, this.local_ind, this.domain)
            cur = cur.next

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_cells(self, SelectorObject selector,
              np.ndarray[np.uint8_t, ndim=1, cast=True] mask):
        cdef int i, j, k, oi
        cdef np.int64_t count = 0
        # pos here is CELL center, not OCT center.
        cdef np.float64_t pos[3]
        cdef int n = mask.shape[0]
        cdef np.float64_t base_dx[3], dx[3]
        cdef int eterm[3]
        assert(self.cont.next == NULL)
        for i in range(3):
            # This is the base_dx, but not the base distance from the center
            # position.  Note that the positions will also all be offset by
            # dx/2.0.
            base_dx[i] = (self.DRE[i] - self.DLE[i])/self.nn[i]
        for oi in range(n):
            if mask[oi] == 0: continue
            o = self.cont.my_octs[oi]
            for i in range(3):
                # This gives the *grid* width for this level
                dx[i] = base_dx[i] / (2 << o.level)
                pos[i] = self.DLE[i] + o.pos[i]*dx[i] + dx[i]/4.0
                dx[i] = dx[i] / 2.0 # This is now the *offset* between cells
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        if o.children[i][j][k] != NULL: continue # child mask
                        count += selector.select_cell(pos, dx, eterm)
                        pos[2] += dx[2]
                    pos[1] += dx[1]
                pos[0] += dx[0]
        return count

cdef class RAMSESOctreeContainer(OctreeContainer):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count(self, np.ndarray[np.uint8_t, ndim=1, cast=True] mask,
                     split = False):
        cdef int n = mask.shape[0]
        cdef int i, dom
        cdef OctAllocationContainer *cur = self.cont
        cdef np.ndarray[np.int64_t, ndim=1] tind
        assert(cur.next == NULL) # Not ready for multiple yet
        cdef np.ndarray[np.int64_t, ndim=1] count
        count = np.zeros(self.max_domain, 'int64')
        for i in range(n):
            if mask[i] == 1:
                count[cur.my_octs[i].domain - 1] += 1
        if split == False:
            return count
        # Now we split them up.
        cdef np.int64_t **arrs = <np.int64_t **> malloc(sizeof(np.int64_t *)
                * self.max_domain)
        cdef np.int64_t *aind = <np.int64_t *> malloc(sizeof(np.int64_t) 
                * self.max_domain)
        tr = []
        for i in range(self.max_domain):
            tind = np.empty(count[i], 'int64')
            tr.append(tind)
            arrs[i] = <np.int64_t *> tind.data
            aind[i] = 0
        for i in range(n):
            if mask[i] == 1:
                dom = cur.my_octs[i].domain - 1
                arrs[dom][aind[dom]] = i
                aind[dom] += 1
        free(arrs)
        free(aind)
        return tr

    def domain_indices(self, np.ndarray[np.uint8_t, ndim=1, cast=True] mask,
                       int domain):
        cdef int n = mask.shape[0]
        cdef int i, p
        p = 0
        cdef OctAllocationContainer *cur = self.cont
        assert(cur.next == NULL) # Not ready for multiple yet
        cdef np.ndarray[np.int64_t, ndim=2] inds
        inds = np.empty((mask.sum(), 2), 'int64')
        for i in range(n):
            if mask[i] == 1:
                inds[p, 0] = cur.my_octs[i].domain
                inds[p, 1] = cur.my_octs[i].ind
                p += 1
        return inds

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def add(self, int curdom, int curlevel, int ng,
            np.ndarray[np.float64_t, ndim=2] pos,
            np.ndarray[np.int64_t, ndim=1] index,
            np.ndarray[np.int64_t, ndim=2] cpumap):
        cdef int level, no, p, i, j, k, ind[3]
        cdef Oct* cur = self.root_mesh[0][0][0]
        cdef np.float64_t pp[3], cp[3], dds[3]
        no = pos.shape[0]
        cdef OctAllocationContainer *cont = self.cont
        # How do we bootstrap ourselves?
        if curdom > self.max_domain: self.max_domain = curdom
        for p in range(no):
            for i in range(3):
                pp[i] = pos[p, i]
                dds[i] = (self.DRE[i] + self.DLE[i])/self.nn[i]
                ind[i] = <np.int64_t> (pp[i]/dds[i])
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
                #print ind[0], ind[1], ind[2]
                next = cur.children[ind[0]][ind[1]][ind[2]]
                if next == NULL:
                    cur.children[ind[0]][ind[1]][ind[2]] = \
                        &cont.my_octs[self.nocts]
                    next = cur.children[ind[0]][ind[1]][ind[2]]
                    next.parent = cur
                    for i in range(3):
                        next.pos[i] = ind[i] + cur.pos[i] * 2
                    next.level = level + 1
                    self.nocts += 1
                cur = next
            cur.domain = curdom
            cur.ind = index[p]
            cur.level = curlevel
