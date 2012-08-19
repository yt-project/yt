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
    if n_cont.my_octs == NULL:
        raise MemoryError
    n_cont.n = n_octs
    n_cont.n_assigned = 0
    for n in range(n_octs):
        oct = &n_cont.my_octs[n]
        oct.parent = NULL
        oct.ind = oct.domain = -1
        oct.local_ind = n + n_cont.offset
        oct.level = -1
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
              np.ndarray[np.uint8_t, ndim=2, cast=True] mask):
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
            o = &cur.my_octs[oi - cur.offset]
            for i in range(8):
                count[o.domain - 1] += mask[oi,i]
        return count

cdef class RAMSESOctreeContainer(OctreeContainer):

    def allocate_domains(self, domain_counts):
        cdef int count, i
        cdef OctAllocationContainer *cur = self.cont
        assert(cur == NULL)
        self.max_domain = len(domain_counts) # 1-indexed
        self.domains = <OctAllocationContainer **> malloc(
            sizeof(OctAllocationContainer *) * len(domain_counts))
        for i, count in enumerate(domain_counts):
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

    def check(self, int curdom):
        cdef int dind, pi
        cdef Oct oct
        cdef OctAllocationContainer *cont = self.domains[curdom - 1]
        cdef int nbad = 0
        for pi in range(cont.n_assigned):
            oct = cont.my_octs[pi]
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        if oct.children[i][j][k] != NULL and \
                           oct.children[i][j][k].level != oct.level + 1:
                            if curdom == 61:
                                print pi, oct.children[i][j][k].level,
                                print oct.level
                            nbad += 1
        print "DOMAIN % 3i HAS % 9i BAD OCTS (%s / %s / %s)" % (curdom, nbad, 
            cont.n - cont.n_assigned, cont.n_assigned, cont.n)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def add(self, int curdom, int curlevel, int ng,
            np.ndarray[np.float64_t, ndim=2] pos,
            int local_domain):
        cdef int level, no, p, i, j, k, ind[3]
        cdef int local = (local_domain == curdom)
        cdef Oct *cur, *next = NULL
        cdef np.float64_t pp[3], cp[3], dds[3]
        no = pos.shape[0] #number of octs
        if curdom > self.max_domain: curdom = local_domain
        cdef OctAllocationContainer *cont = self.domains[curdom - 1]
        cdef int initial = cont.n_assigned
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
                if curlevel != 0:
                    raise RuntimeError
                cur = &cont.my_octs[cont.n_assigned]
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
                next = cur.children[ind[0]][ind[1]][ind[2]]
                if next == NULL:
                    next = &cont.my_octs[cont.n_assigned]
                    cur.children[ind[0]][ind[1]][ind[2]] = next
                    cont.n_assigned += 1
                    next.parent = cur
                    for i in range(3):
                        next.pos[i] = ind[i] + (cur.pos[i] << 1)
                    next.level = level + 1
                    self.nocts += 1
                cur = next
            cur.domain = curdom
            if local == 1:
                cur.ind = p
            cur.level = curlevel
        return cont.n_assigned - initial

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def icoords(self, int domain_id,
                np.ndarray[np.uint8_t, ndim=2, cast=True] mask,
                np.int64_t cell_count,
                np.ndarray[np.int64_t, ndim=1] level_counts):
        # Wham, bam, it's a scam
        cdef np.int64_t i, j, k, oi, ci, n, ii, level
        cdef OctAllocationContainer *cur = self.domains[domain_id - 1]
        cdef Oct *o
        n = mask.shape[0]
        cdef np.ndarray[np.int64_t, ndim=2] coords
        coords = np.empty((cell_count, 3), dtype="int64")
        for oi in range(cur.n):
            o = &cur.my_octs[oi]
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        ii = ((k*2)+j)*2+i
                        if mask[o.local_ind, ii] == 0: continue
                        ci = level_counts[o.level]
                        coords[ci, 0] = (o.pos[0] << 1) + i
                        coords[ci, 1] = (o.pos[1] << 1) + j
                        coords[ci, 2] = (o.pos[2] << 1) + k
                        level_counts[o.level] += 1
        return coords

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def ires(self, int domain_id,
                np.ndarray[np.uint8_t, ndim=2, cast=True] mask,
                np.int64_t cell_count,
                np.ndarray[np.int64_t, ndim=1] level_counts):
        # Wham, bam, it's a scam
        cdef np.int64_t i, j, k, oi, ci, n
        cdef OctAllocationContainer *cur = self.domains[domain_id - 1]
        cdef Oct *o
        n = mask.shape[0]
        cdef np.ndarray[np.int64_t, ndim=1] levels
        levels = np.empty(cell_count, dtype="int64")
        ci = 0
        for oi in range(cur.n):
            o = &cur.my_octs[oi]
            for i in range(8):
                if mask[oi + cur.offset, i] == 0: continue
                ci = level_counts[o.level]
                levels[ci] = o.level
                level_counts[o.level] += 1
        return levels

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_levels(self, int max_level, int domain_id,
                     np.ndarray[np.uint8_t, ndim=2, cast=True] mask):
        cdef np.ndarray[np.int64_t, ndim=1] level_count
        cdef OctAllocationContainer *cur = self.domains[domain_id - 1]
        cdef Oct *o
        cdef int oi, i
        level_count = np.zeros(max_level+1, 'int64')
        for oi in range(cur.n):
            o = &cur.my_octs[oi]
            for i in range(8):
                if mask[o.local_ind, i] == 0: continue
                level_count[o.level] += 1
        return level_count

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fcoords(self, int domain_id,
                np.ndarray[np.uint8_t, ndim=2, cast=True] mask,
                np.int64_t cell_count,
                np.ndarray[np.int64_t, ndim=1] level_counts):
        # Wham, bam, it's a scam
        cdef np.int64_t i, j, k, oi, ci, n, ii
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
            # dx/2.0.  This is also for *oct grids*, not cells.
            base_dx[i] = (self.DRE[i] - self.DLE[i])/self.nn[i]
        for oi in range(cur.n):
            o = &cur.my_octs[oi]
            for i in range(3):
                # This gives the *grid* width for this level
                dx[i] = base_dx[i] / (1 << o.level)
                # o.pos is the *grid* index, so pos[i] is the center of the
                # first cell in the grid
                pos[i] = self.DLE[i] + o.pos[i]*dx[i] + dx[i]/4.0
                dx[i] = dx[i] / 2.0 # This is now the *offset* 
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        ii = ((k*2)+j)*2+i
                        if mask[o.local_ind, ii] == 0: continue
                        ci = level_counts[o.level]
                        coords[ci, 0] = pos[0] + dx[0] * i
                        coords[ci, 1] = pos[1] + dx[1] * j
                        coords[ci, 2] = pos[2] + dx[2] * k
                        level_counts[o.level] += 1
        return coords

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_level(self, int domain, int level, dest_fields, source_fields,
                   np.ndarray[np.uint8_t, ndim=2, cast=True] mask, int offset):
        cdef np.ndarray[np.float64_t, ndim=2] source
        cdef np.ndarray[np.float64_t, ndim=1] dest
        cdef OctAllocationContainer *dom = self.domains[domain - 1]
        cdef Oct *o
        cdef int n
        cdef int i, j, k, ii
        cdef int local_pos, local_filled
        cdef np.float64_t val
        for key in dest_fields:
            local_filled = 0
            dest = dest_fields[key]
            source = source_fields[key]
            for n in range(dom.n):
                o = &dom.my_octs[n]
                if o.level != level: continue
                for i in range(2):
                    for j in range(2):
                        for k in range(2):
                            ii = ((k*2)+j)*2+i
                            if mask[o.local_ind, ii] == 0: continue
                            dest[local_filled + offset] = source[o.ind, ii]
                            local_filled += 1
        return local_filled

cdef class ParticleOctreeContainer(OctreeContainer):
    cdef ParticleArrays *first_sd
    cdef ParticleArrays *last_sd
    cdef Oct** oct_list

    def __dealloc__(self):
        cdef i, j, k
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                for k in range(self.nn[2]):
                    self.visit_free(self.root_mesh[i][j][k])

    cdef void visit_free(self, Oct *o):
        cdef int i, j, k
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if o.children[i][j][k] == NULL: continue
                    self.visit_free(o.children[i][j][k])
        if o.sd.np >= 0:
            if o.sd.pos != NULL:
                for i in range(3):
                    free(o.sd.pos[i])
                free(o.sd.pos)
            free(o.sd.domain_id)
        free(o)

    def allocate_domains(self, domain_counts):
        cdef int count, i

    def finalize(self):
        self.oct_list = <Oct**> malloc(sizeof(Oct*)*self.nocts)
        cdef i = 0
        cdef ParticleArrays *c = self.first_sd
        while c != NULL:
            self.oct_list[i] = c.oct
            if c.np >= 0:
                for j in range(3):
                    free(c.pos[j])
                free(c.pos)
                c.pos = NULL
                # We should also include a shortening of the domain IDs here
            c = c.next
            i += 1

    cdef Oct* allocate_oct(self):
        self.nocts += 1
        cdef Oct *my_oct = <Oct*> malloc(sizeof(Oct))
        cdef ParticleArrays *sd = <ParticleArrays*> \
            malloc(sizeof(ParticleArrays))
        cdef int i, j, k
        my_oct.ind = my_oct.domain = -1
        my_oct.local_ind = self.nocts - 1
        my_oct.pos[0] = my_oct.pos[1] = my_oct.pos[2] = -1
        my_oct.level = -1
        my_oct.sd = sd
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    my_oct.children[i][j][k] = NULL
        my_oct.parent = NULL
        if self.first_sd == NULL:
            self.first_sd = sd
        if self.last_sd != NULL:
            self.last_sd.next = sd
        self.last_sd = sd
        sd.oct = my_oct
        sd.next = NULL
        sd.domain_id = <np.int64_t *> malloc(sizeof(np.int64_t) * 32)
        sd.pos = <np.float64_t **> malloc(sizeof(np.float64_t*) * 3)
        for i in range(3):
            sd.pos[i] = <np.float64_t *> malloc(sizeof(np.float64_t) * 32)
        for i in range(32):
            sd.pos[0][i] = sd.pos[1][i] = sd.pos[2][i] = 0.0
            sd.domain_id[i] = -1
        sd.np = 0
        return my_oct

    def linearly_count(self):
        cdef np.int64_t total = 0
        cdef ParticleArrays *c = self.first_sd
        while c != NULL:
            total += 1
            c = c.next
        return total

    def add(self, np.ndarray[np.float64_t, ndim=2] pos, np.int64_t domain_id):
        cdef int no = pos.shape[0]
        cdef int p, i, level
        cdef np.float64_t dds[3], cp[3], pp[3]
        cdef int ind[3]
        self.max_domain = max(self.max_domain, domain_id)
        for p in range(no):
            level = 0
            for i in range(3):
                pp[i] = pos[p, i]
                dds[i] = (self.DRE[i] + self.DLE[i])/self.nn[i]
                ind[i] = <np.int64_t> (pp[i]/dds[i])
                cp[i] = (ind[i] + 0.5) * dds[i]
            cur = self.root_mesh[ind[0]][ind[1]][ind[2]]
            if cur == NULL:
                cur = self.allocate_oct()
                self.root_mesh[ind[0]][ind[1]][ind[2]] = cur
                for i in range(3):
                    cur.pos[i] = ind[i] # root level
            if cur.sd.np == 32:
                self.refine_oct(cur, cp)
            while cur.sd.np < 0:
                for i in range(3):
                    dds[i] = dds[i] / 2.0
                    if cp[i] > pp[i]:
                        ind[i] = 0
                        cp[i] -= dds[i] / 2.0
                    else:
                        ind[i] = 1
                        cp[i] += dds[i]/2.0
                cur = cur.children[ind[0]][ind[1]][ind[2]]
                level += 1
                if cur.sd.np == 32:
                    self.refine_oct(cur, cp)
            # Now we copy in our particle 
            pi = cur.sd.np
            for i in range(3):
                cur.sd.pos[i][pi] = pp[i]
            cur.sd.domain_id[pi] = domain_id
            cur.sd.np += 1

    cdef void refine_oct(self, Oct *o, np.float64_t pos[3]):
        cdef int i, j, k, m, ind[3]
        cdef Oct *noct
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    noct = self.allocate_oct()
                    noct.pos[0] = o.pos[0] << 1 + i
                    noct.pos[1] = o.pos[1] << 1 + j
                    noct.pos[2] = o.pos[2] << 1 + k
                    o.children[i][j][k] = noct
        for m in range(32):
            for i in range(3):
                if o.sd.pos[i][m] < pos[i]:
                    ind[i] = 0
                else:
                    ind[i] = 1
            noct = o.children[ind[0]][ind[1]][ind[2]]
            k = noct.sd.np
            for i in range(3):
                noct.sd.pos[i][k] = o.sd.pos[i][m]
            noct.sd.domain_id[k] = o.sd.domain_id[k]
            noct.sd.np += 1
        o.sd.np = -1
        for i in range(3):
            free(o.sd.pos[i])
        free(o.sd.domain_id)
        free(o.sd.pos)

    def recursively_count(self):
        cdef int i, j, k
        cdef np.int64_t counts[128]
        for i in range(128): counts[i] = 0
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                for k in range(self.nn[2]):
                    if self.root_mesh[i][j][k] != NULL:
                        self.visit(self.root_mesh[i][j][k], counts)
        level_counts = {}
        for i in range(128):
            if counts[i] == 0: break
            level_counts[i] = counts[i]
        return level_counts
        
    cdef visit(self, Oct *o, np.int64_t *counts, level = 0):
        cdef int i, j, k
        counts[level] += 1
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if o.children[i][j][k] != NULL:
                        self.visit(o.children[i][j][k], counts, level + 1)
        return

    def domain_identify(self, np.ndarray[np.uint8_t, ndim=2, cast=True] mask):
        cdef int i, oi, m
        cdef Oct *o
        cdef np.ndarray[np.uint8_t, ndim=1, cast=True] dmask
        dmask = np.zeros(self.max_domain+1, dtype='uint8')
        for oi in range(self.nocts):
            m = 0
            o = self.oct_list[oi]
            if o.sd.np <= 0: continue
            for i in range(8):
                if mask[oi, i] == 1:
                    m = 1
                    break
            if m == 0: continue
            for i in range(o.sd.np):
                dmask[o.sd.domain_id[i]] = 1
        return dmask.astype("bool")
