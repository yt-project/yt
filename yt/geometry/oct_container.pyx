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
    if prev == NULL:
        n_cont.offset = 0
    else:
        n_cont.offset = prev.offset + prev.n
    n_cont.my_octs = <Oct *> malloc(sizeof(Oct) * n_octs)
    n_cont.n = n_octs
    n_cont.n_assigned = 0
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

    def __init__(self, domain_dimensions, domain_left_edge, domain_right_edge):
        self.nn[0], self.nn[1], self.nn[2] = domain_dimensions
        cdef int i, j, k, p
        self.max_domain = -1
        p = 0
        self.nocts = 0 # Increment when initialized
        self.root_mesh = <Oct****> malloc(sizeof(void*) * self.nn[0])
        for i in range(3):
            self.nn[i] = domain_dimensions[i]
        for i in range(self.nn[0]):
            self.root_mesh[i] = <Oct ***> malloc(sizeof(void*) * self.nn[1])
            for j in range(self.nn[1]):
                self.root_mesh[i][j] = <Oct **> malloc(sizeof(void*) * self.nn[2])
                for k in range(self.nn[2]):
                    self.root_mesh[i][j][k] = NULL
        # We don't initialize the octs yet
        for i in range(3):
            self.DLE[i] = domain_left_edge[i]
            self.DRE[i] = domain_right_edge[i]

    def __dealloc__(self):
        free_octs(self.cont)
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                if self.root_mesh[i][j] == NULL: continue
                free(self.root_mesh[i][j])
            if self.root_mesh[i] == NULL: continue
            free(self.root_mesh[i])
        if self.root_mesh == NULL: return
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
        # pos here is CELL center, not OCT center.
        cdef np.float64_t pos[3]
        cdef int n = mask.shape[0]
        cdef np.float64_t base_dx[3], dx[3]
        cdef int eterm[3]
        cdef np.ndarray[np.int64_t, ndim=1] count
        count = np.zeros(self.max_domain, 'int64')
        for i in range(3):
            # This is the base_dx, but not the base distance from the center
            # position.  Note that the positions will also all be offset by
            # dx/2.0.
            base_dx[i] = (self.DRE[i] - self.DLE[i])/self.nn[i]
        cur = self.cont
        for oi in range(n):
            if oi - cur.offset >= cur.n:
                cur = cur.next
            if mask[oi] == 0: continue
            o = &cur.my_octs[oi - cur.offset]
            for i in range(3):
                # This gives the *grid* width for this level
                dx[i] = base_dx[i] / (2 << o.level)
                pos[i] = self.DLE[i] + o.pos[i]*dx[i] + dx[i]/4.0
                dx[i] = dx[i] / 2.0 # This is now the *offset* 
            for i in range(2):
                pos[1] = self.DLE[1] + o.pos[1]*dx[1] + dx[1]/4.0
                for j in range(2):
                    pos[2] = self.DLE[2] + o.pos[2]*dx[2] + dx[2]/4.0
                    for k in range(2):
                        if o.children[i][j][k] == NULL: 
                            count[o.domain - 1] += \
                                selector.select_cell(pos, dx, eterm)
                        pos[2] += dx[2]
                    pos[1] += dx[1]
                pos[0] += dx[0]
        return count

cdef class RAMSESOctreeContainer(OctreeContainer):

    def allocate_domains(self, domain_counts):
        cdef int count, i
        cdef OctAllocationContainer *cur = self.cont
        assert(cur == NULL)
        self.max_domain = len(domain_counts) # 1-indexed
        self.domains = <OctAllocationContainer **> malloc(
            sizeof(OctAllocationContainer *) * len(domain_counts))
        for i,count in enumerate(domain_counts):
            cur = allocate_octs(count, cur)
            if self.cont == NULL: self.cont = cur
            self.domains[i] = cur
        
    def __dealloc__(self):
        # This gets called BEFORE the superclass deallocation.  But, both get
        # called.
        if self.domains != NULL: free(self.domains)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count(self, np.ndarray[np.uint8_t, ndim=1, cast=True] mask,
                     split = False):
        cdef int n = mask.shape[0]
        cdef int i, dom
        cdef OctAllocationContainer *cur
        cdef np.ndarray[np.int64_t, ndim=1] count
        count = np.zeros(self.max_domain, 'int64')
        # This is the idiom for iterating over many containers.
        cur = self.cont
        for i in range(n):
            if i - cur.offset >= cur.n: cur = cur.next
            if mask[i] == 1:
                count[cur.my_octs[i - cur.offset].domain - 1] += 1
        return count

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def add(self, int curdom, int curlevel, int ng,
            np.ndarray[np.float64_t, ndim=2] pos,
            np.ndarray[np.int64_t, ndim=1] findices,
            np.ndarray[np.int64_t, ndim=2] cpumap):
        cdef int level, no, p, i, j, k, ind[3]
        cdef Oct* cur = self.root_mesh[0][0][0]
        cdef np.float64_t pp[3], cp[3], dds[3]
        no = pos.shape[0] #number of octs
        cdef OctAllocationContainer *cont = self.domains[curdom - 1]
        # How do we bootstrap ourselves?
        for p in range(no):
            #for every oct we're trying to add find the 
            #floating point unitary position on this level
            for i in range(3):
                pp[i] = pos[p, i]
                dds[i] = (self.DRE[i] + self.DLE[i])/self.nn[i]
                ind[i] = <np.int64_t> (pp[i]/dds[i])
                cp[i] = (ind[i] + 0.5) * dds[i]
            cur = self.root_mesh[ind[0]][ind[1]][ind[2]]
            if cur == NULL:
                if curlevel != 0: raise RuntimeError
                cur = &cont.my_octs[cont.n_assigned]
                cur.local_ind = cont.offset + cont.n_assigned
                cur.parent = NULL
                cur.level = 0
                for i in range(3):
                    cur.pos[i] = ind[i]
                cont.n_assigned += 1
                self.nocts += 1
                self.root_mesh[ind[0]][ind[1]][ind[2]] = cur
            # Now we find the location we want
            # Note that RAMSES I think 1-findiceses levels, but we don't.
            for level in range(curlevel):
                # At every level, find the cell this oct
                # lives inside
                for i in range(3):
                    #as we get deeper, oct size halves
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
                        &cont.my_octs[cont.n_assigned]
                    next = cur.children[ind[0]][ind[1]][ind[2]]
                    next.local_ind = cont.offset + cont.n_assigned
                    next.parent = cur
                    for i in range(3):
                        next.pos[i] = ind[i] + (cur.pos[i] << 1)
                    next.level = level + 1
                    cont.n_assigned += 1
                    self.nocts += 1
                cur = next
            cur.domain = curdom
            cur.ind = findices[p]
            cur.level = curlevel

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def icoords(self, int domain_id,
                np.ndarray[np.uint8_t, ndim=1, cast=True] mask,
                np.int64_t cell_count):
        # Wham, bam, it's a scam
        cdef np.int64_t i, j, k, oi, ci, n
        cdef OctAllocationContainer *cur = self.domains[domain_id - 1]
        cdef Oct *o
        n = mask.shape[0]
        cdef np.ndarray[np.int64_t, ndim=2] coords
        coords = np.empty((cell_count, 3), dtype="int64")

        ci = 0
        for oi in range(cur.n):
            if mask[oi + cur.offset] == 0: continue
            o = &cur.my_octs[oi]
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        if o.children[i][j][k] != NULL: continue
                        coords[ci, 0] = (o.pos[0] << 1) + i
                        coords[ci, 1] = (o.pos[1] << 1) + j
                        coords[ci, 2] = (o.pos[2] << 1) + k
                        ci += 1
        return coords

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def ires(self, int domain_id,
                np.ndarray[np.uint8_t, ndim=1, cast=True] mask,
                np.int64_t cell_count):
        # Wham, bam, it's a scam
        cdef np.int64_t i, j, k, oi, ci, n
        cdef OctAllocationContainer *cur = self.domains[domain_id - 1]
        cdef Oct *o
        n = mask.shape[0]
        cdef np.ndarray[np.int64_t, ndim=1] levels
        levels = np.empty(cell_count, dtype="int64")
        ci = 0
        for oi in range(cur.n):
            if mask[oi + cur.offset] == 0: continue
            o = &cur.my_octs[oi]
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        if o.children[i][j][k] != NULL: continue
                        levels[ci] = o.level
                        ci += 1
        return levels

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fcoords(self, int domain_id,
                np.ndarray[np.uint8_t, ndim=1, cast=True] mask,
                np.int64_t cell_count):
        # Wham, bam, it's a scam
        cdef np.int64_t i, j, k, oi, ci, n
        cdef OctAllocationContainer *cur = self.domains[domain_id - 1]
        cdef Oct *o
        cdef np.float64_t pos[3]
        cdef np.float64_t base_dx[3], dx[3]
        n = mask.shape[0]
        cdef np.ndarray[np.float64_t, ndim=2] coords
        coords = np.empty((cell_count, 3), dtype="float64")
        for i in range(3):
            # This is the base_dx, but not the base distance from the center
            # position.  Note that the positions will also all be offset by
            # dx/2.0.
            base_dx[i] = (self.DRE[i] - self.DLE[i])/self.nn[i]
        ci = 0
        for oi in range(cur.n):
            if mask[oi + cur.offset] == 0: continue
            o = &cur.my_octs[oi]
            for i in range(3):
                # This gives the *grid* width for this level
                dx[i] = base_dx[i] / (2 << o.level)
                pos[i] = self.DLE[i] + o.pos[i]*dx[i] + dx[i]/4.0
                dx[i] = dx[i] / 2.0 # This is now the *offset* 
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        if o.children[i][j][k] != NULL: continue
                        coords[ci, 0] = pos[0] + dx[0] * i
                        coords[ci, 1] = pos[1] + dx[1] * j
                        coords[ci, 2] = pos[2] + dx[2] * k
                        ci += 1
        return coords
