"""
Oct container

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Author: Christopher Moody <chris.e.moody@gmail.com>
Affiliation: UC Santa Cruz
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

from libc.stdlib cimport malloc, free, qsort
from libc.math cimport floor
cimport numpy as np
import numpy as np
from oct_container cimport Oct, OctAllocationContainer, \
    OctreeContainer, ORDER_MAX
from selection_routines cimport SelectorObject, \
    OctVisitorData, oct_visitor_function
cimport cython

ORDER_MAX = 20
_ORDER_MAX = ORDER_MAX

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
        oct.file_ind = oct.domain = -1
        oct.domain_ind = n + n_cont.offset
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
#   * Note that this does not allow OctAllocationContainer to exactly be a
#     chunk, but it is close.  For IO chunking, we can theoretically examine
#     those octs that live inside a given allocator.

cdef class OctreeContainer:

    def __init__(self, oct_domain_dimensions, domain_left_edge, domain_right_edge):
        # This will just initialize the root mesh octs
        cdef int i, j, k, p
        for i in range(3):
            self.nn[i] = oct_domain_dimensions[i]
        self.max_domain = -1
        p = 0
        self.nocts = 0 # Increment when initialized
        self.root_mesh = <Oct****> malloc(sizeof(void*) * self.nn[0])
        for i in range(self.nn[0]):
            self.root_mesh[i] = <Oct ***> malloc(sizeof(void*) * self.nn[1])
            for j in range(self.nn[1]):
                self.root_mesh[i][j] = <Oct **> malloc(sizeof(void*) * self.nn[2])
                for k in range(self.nn[2]):
                    self.root_mesh[i][j][k] = NULL
        # We don't initialize the octs yet
        for i in range(3):
            self.DLE[i] = domain_left_edge[i] #0
            self.DRE[i] = domain_right_edge[i] #num_grid

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
        #Get the next oct, will traverse domains
        #Note that oct containers can be sorted 
        #so that consecutive octs are on the same domain
        cdef OctAllocationContainer *cur = self.cont
        cdef Oct *this
        cdef int i
        while cur != NULL:
            for i in range(cur.n_assigned):
                this = &cur.my_octs[i]
                yield (this.file_ind, this.domain_ind, this.domain)
            cur = cur.next

    cdef void visit_all_octs(self, SelectorObject selector,
                        oct_visitor_function *func,
                        OctVisitorData *data):
        cdef int i, j, k, n
        cdef np.float64_t pos[3], dds[3]
        # This dds is the oct-width
        for i in range(3):
            dds[i] = (self.DRE[i] - self.DLE[i]) / self.nn[i]
        # Pos is the center of the octs
        pos[0] = self.DLE[0] + dds[0]/2.0
        for i in range(self.nn[0]):
            pos[1] = self.DLE[1] + dds[1]/2.0
            for j in range(self.nn[1]):
                pos[2] = self.DLE[2] + dds[2]/2.0
                for k in range(self.nn[2]):
                    if self.root_mesh[i][j][k] == NULL: continue
                    selector.recursively_visit_octs(
                        self.root_mesh[i][j][k],
                        pos, dds, 0, func, data)
                    pos[2] += dds[2]
                pos[1] += dds[1]
            pos[0] += dds[0]

    cdef void oct_bounds(self, Oct *o, np.float64_t *corner, np.float64_t *size):
        cdef int i
        for i in range(3):
            size[i] = (self.DRE[i] - self.DLE[i]) / (self.nn[i] << o.level)
            corner[i] = o.pos[i] * size[i] + self.DLE[i]

    cdef np.int64_t get_domain_offset(self, int domain_id):
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Oct *get(self, np.float64_t ppos[3], OctInfo *oinfo = NULL):
        #Given a floating point position, retrieve the most
        #refined oct at that time
        cdef np.int64_t ind[3]
        cdef np.float64_t dds[3], cp[3], pp[3]
        cdef Oct *cur
        cdef int i
        for i in range(3):
            dds[i] = (self.DRE[i] - self.DLE[i])/self.nn[i]
            ind[i] = <np.int64_t> ((ppos[i] - self.DLE[i])/dds[i])
            cp[i] = (ind[i] + 0.5) * dds[i] + self.DLE[i]
        next = self.root_mesh[ind[0]][ind[1]][ind[2]]
        # We want to stop recursing when there's nowhere else to go
        while next != NULL:
            cur = next
            for i in range(3):
                dds[i] = dds[i] / 2.0
                if cp[i] > ppos[i]:
                    ind[i] = 0
                    cp[i] -= dds[i] / 2.0
                else:
                    ind[i] = 1
                    cp[i] += dds[i]/2.0
            next = cur.children[ind[0]][ind[1]][ind[2]]
        if oinfo == NULL: return cur
        for i in range(3):
            # This will happen *after* we quit out, so we need to back out the
            # last change to cp
            if ind[i] == 1:
                cp[i] -= dds[i]/2.0 # Now centered
            else:
                cp[i] += dds[i]/2.0
            # We don't need to change dds[i] as it has been halved from the
            # oct width, thus making it already the cell width
            oinfo.dds[i] = dds[i] # Cell width
            oinfo.left_edge[i] = cp[i] - dds[i] # Center minus dds
        return cur

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_cells(self, SelectorObject selector,
              np.ndarray[np.uint8_t, ndim=2, cast=True] mask):
        cdef int i, j, k
        cdef np.int64_t oi
        # pos here is CELL center, not OCT center.
        cdef np.float64_t pos[3]
        cdef int n = mask.shape[0]
        cdef np.ndarray[np.int64_t, ndim=1] count
        count = np.zeros(self.max_domain, 'int64')
        # 
        cur = self.cont
        for oi in range(n):
            if oi - cur.offset >= cur.n_assigned:
                cur = cur.next
            o = &cur.my_octs[oi - cur.offset]
            for i in range(8):
                count[o.domain - 1] += mask[o.domain_ind,i]
        return count

    @cython.boundscheck(True)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_leaves(self, np.ndarray[np.uint8_t, ndim=2, cast=True] mask):
        # Modified to work when not all octs are assigned
        cdef int i, j, k, ii
        cdef np.int64_t oi
        # pos here is CELL center, not OCT center.
        cdef np.float64_t pos[3]
        cdef int n = mask.shape[0]
        cdef np.ndarray[np.int64_t, ndim=1] count
        count = np.zeros(self.max_domain, 'int64')
        # 
        cur = self.cont
        for oi in range(n):
            if oi - cur.offset >= cur.n_assigned:
                cur = cur.next
                if cur == NULL:
                    break
            o = &cur.my_octs[oi - cur.offset]
            # skip if unassigned
            if o == NULL:
                continue
            if o.domain == -1: 
                continue
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        if o.children[i][j][k] == NULL:
                            ii = ((k*2)+j)*2+i
                            count[o.domain - 1] += mask[o.domain_ind,ii]
        return count

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void neighbors(self, Oct* o, Oct* neighbors[27]):
        #Get 3x3x3 neighbors, although the 1,1,1 oct is the
        #central one. 
        #Return an array of Octs
        cdef np.int64_t curopos[3]
        cdef np.int64_t curnpos[3]
        cdef np.int64_t npos[3]
        cdef int i, j, k, ni, nj, nk, ind[3], nn, dl, skip
        cdef np.float64_t dds[3], cp[3], pp[3]
        cdef Oct* candidate
        for i in range(27): neighbors[i] = NULL
        nn = 0
        for ni in range(3):
            for nj in range(3):
                for nk in range(3):
                    if ni == nj == nk == 1:
                        neighbors[nn] = o
                        nn += 1
                        continue
                    npos[0] = o.pos[0] + (ni - 1)
                    npos[1] = o.pos[1] + (nj - 1)
                    npos[2] = o.pos[2] + (nk - 1)
                    for i in range(3):
                        # Periodicity
                        if npos[i] == -1:
                            npos[i] = (self.nn[i]  << o.level) - 1
                        elif npos[i] == (self.nn[i] << o.level):
                            npos[i] = 0
                        curopos[i] = o.pos[i]
                        curnpos[i] = npos[i] 
                    # Now we have our neighbor position and a safe place to
                    # keep it.  curnpos will be the root index of the neighbor
                    # at a given level, and npos will be constant.  curopos is
                    # the candidate root at a level.
                    candidate = o
                    while candidate != NULL:
                        if ((curopos[0] == curnpos[0]) and 
                            (curopos[1] == curnpos[1]) and
                            (curopos[2] == curnpos[2])):
                            break
                        # This one doesn't meet it, so we pop up a level.
                        # First we update our positions, then we update our
                        # candidate.
                        for i in range(3):
                            # We strip a digit off the right
                            curopos[i] = (curopos[i] >> 1)
                            curnpos[i] = (curnpos[i] >> 1)
                        # Now we update to the candidate's parent, which should
                        # have a matching position to curopos[]
                        candidate = candidate.parent
                    if candidate == NULL:
                        # Worst case scenario
                        for i in range(3):
                            ind[i] = (npos[i] >> (o.level))
                        candidate = self.root_mesh[ind[0]][ind[1]][ind[2]]
                    # Now we have the common root, which may be NULL
                    while candidate.level < o.level:
                        dl = o.level - (candidate.level + 1)
                        for i in range(3):
                            ind[i] = (npos[i] >> dl) & 1
                        if candidate.children[0][0][0] == NULL: break
                        candidate = candidate.children[ind[0]][ind[1]][ind[2]]
                    neighbors[nn] = candidate
                    nn += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def get_neighbor_boundaries(self, oppos):
        cdef int i, ii
        cdef np.float64_t ppos[3]
        for i in range(3):
            ppos[i] = oppos[i]
        cdef Oct *main = self.get(ppos)
        cdef Oct* neighbors[27]
        self.neighbors(main, neighbors)
        cdef np.ndarray[np.float64_t, ndim=2] bounds
        cdef np.float64_t corner[3], size[3]
        bounds = np.zeros((27,6), dtype="float64")
        tnp = 0
        for i in range(27):
            self.oct_bounds(neighbors[i], corner, size)
            for ii in range(3):
                bounds[i, ii] = corner[ii]
                bounds[i, 3+ii] = size[ii]
        return bounds

cdef class RAMSESOctreeContainer(OctreeContainer):

    cdef np.int64_t get_domain_offset(self, int domain_id):
        cdef OctAllocationContainer *cont = self.domains[domain_id - 1]
        return cont.offset

    cdef Oct* next_root(self, int domain_id, int ind[3]):
        cdef Oct *next = self.root_mesh[ind[0]][ind[1]][ind[2]]
        if next != NULL: return next
        cdef OctAllocationContainer *cont = self.domains[domain_id - 1]
        if cont.n_assigned >= cont.n: raise RuntimeError
        next = &cont.my_octs[cont.n_assigned]
        cont.n_assigned += 1
        self.root_mesh[ind[0]][ind[1]][ind[2]] = next
        next.parent = NULL
        next.level = 0
        for i in range(3):
            next.pos[i] = ind[i]
        self.nocts += 1
        return next

    cdef Oct* next_child(self, int domain_id, int ind[3], Oct *parent):
        cdef Oct *next = parent.children[ind[0]][ind[1]][ind[2]]
        if next != NULL: return next
        cdef OctAllocationContainer *cont = self.domains[domain_id - 1]
        if cont.n_assigned >= cont.n: raise RuntimeError
        next = &cont.my_octs[cont.n_assigned]
        cont.n_assigned += 1
        parent.children[ind[0]][ind[1]][ind[2]] = next
        next.parent = parent
        next.level = parent.level + 1
        for i in range(3):
            next.pos[i] = ind[i] + (parent.pos[i] << 1)
        self.nocts += 1
        return next

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
            if i - cur.offset >= cur.n_assigned: cur = cur.next
            if mask[i] == 1:
                count[cur.my_octs[i - cur.offset].domain - 1] += 1
        return count

    def domain_and(self, np.ndarray[np.uint8_t, ndim=2, cast=True] mask,
                   int domain_id):
        cdef np.int64_t i, oi, n,  use
        cdef OctAllocationContainer *cur = self.domains[domain_id - 1]
        cdef Oct *o
        cdef np.ndarray[np.uint8_t, ndim=2] m2 = \
                np.zeros((mask.shape[0], 8), 'uint8')
        n = mask.shape[0]
        for oi in range(cur.n_assigned):
            o = &cur.my_octs[oi]
            use = 0
            for i in range(8):
                m2[o.domain_ind, i] = mask[o.domain_ind, i]
        return m2 # NOTE: This is uint8_t

    def domain_mask(self,
                    # mask is the base selector's *global* mask
                    np.ndarray[np.uint8_t, ndim=2, cast=True] mask,
                    int domain_id):
        # What distinguishes this one from domain_and is that we have a mask,
        # which covers the whole domain, but our output will only be of a much
        # smaller subset of octs that belong to a given domain *and* the mask.
        # Note also that typically when something calls domain_and, they will 
        # use a logical_any along the oct axis.  Here we don't do that.
        # Note also that we change the shape of the returned array.
        cdef np.int64_t i, j, k, oi, n, nm, use
        cdef OctAllocationContainer *cur = self.domains[domain_id - 1]
        cdef Oct *o
        n = mask.shape[0]
        nm = 0
        for oi in range(cur.n_assigned):
            o = &cur.my_octs[oi]
            use = 0
            for i in range(8):
                if mask[o.domain_ind, i] == 1: use = 1
            nm += use
        cdef np.ndarray[np.uint8_t, ndim=4] m2 = \
                np.zeros((2, 2, 2, nm), 'uint8')
        nm = 0
        for oi in range(cur.n_assigned):
            o = &cur.my_octs[oi]
            use = 0
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        ii = ((k*2)+j)*2+i
                        if mask[o.domain_ind, ii] == 0: continue
                        use = m2[i, j, k, nm] = 1
            nm += use
        return m2.astype("bool")

    def domain_ind(self,
                    # mask is the base selector's *global* mask
                    np.ndarray[np.uint8_t, ndim=2, cast=True] mask,
                    int domain_id):
        cdef np.int64_t i, j, k, oi, noct, n, nm, use, offset
        cdef OctAllocationContainer *cur = self.domains[domain_id - 1]
        cdef Oct *o
        cdef np.ndarray[np.int64_t, ndim=1] ind = np.zeros(cur.n, 'int64') - 1
        nm = 0
        for oi in range(cur.n):
            o = &cur.my_octs[oi]
            use = 0
            for i in range(8):
                if mask[o.domain_ind, i] == 1: use = 1
            if use == 1:
                ind[o.domain_ind - cur.offset] = nm
            nm += use
        return ind

    def check(self, int curdom, int print_all = 0):
        cdef int dind, pi
        cdef Oct oct
        cdef OctAllocationContainer *cont = self.domains[curdom - 1]
        cdef int nbad = 0
        cdef int nmissed = 0
        cdef int unassigned = 0
        for pi in range(cont.n_assigned):
            oct = cont.my_octs[pi]
            if print_all==1:
                print pi, oct.level, oct.domain,
                print oct.pos[0],oct.pos[1],oct.pos[2]
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        if oct.children[i][j][k] != NULL and \
                           oct.children[i][j][k].level != oct.level + 1:
                            nbad += 1
                        if oct.domain != curdom:
                            print curdom, oct.domain
                            nmissed += 1
                        if oct.domain == -1:
                            unassigned += 1
        print "DOMAIN % 3i HAS % 9i BAD OCTS (%s / %s / %s)" % (curdom, nbad, 
            cont.n - cont.n_assigned, cont.n_assigned, cont.n)
        print "DOMAIN % 3i HAS % 9i MISSED OCTS" % (curdom, nmissed)
        print "DOMAIN % 3i HAS % 9i UNASSIGNED OCTS" % (curdom, unassigned)

    def check_refinement(self, int curdom):
        cdef int pi, i, j, k, some_refined, some_unrefined
        cdef Oct *oct
        cdef int bad = 0
        cdef OctAllocationContainer *cont = self.domains[curdom - 1]
        for pi in range(cont.n_assigned):
            oct = &cont.my_octs[pi]
            some_unrefined = 0
            some_refined = 0
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        if oct.children[i][j][k] == NULL:
                            some_unrefined = 1
                        else:
                            some_refined = 1
            if some_unrefined == some_refined == 1:
                #print "BAD", oct.file_ind, oct.domain_ind
                bad += 1
                if curdom == 10 or curdom == 72:
                    for i in range(2):
                        for j in range(2):
                            for k in range(2):
                                print (oct.children[i][j][k] == NULL),
                    print
        print "BAD TOTAL", curdom, bad, cont.n_assigned

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def add(self, int curdom, int curlevel, int ng,
            np.ndarray[np.float64_t, ndim=2] pos,
            int local_domain, int skip_boundary = 1):
        cdef int level, no, p, i, j, k, ind[3]
        cdef int local = (local_domain == curdom)
        cdef Oct *cur, *next = NULL
        cdef np.float64_t pp[3], cp[3], dds[3]
        no = pos.shape[0] #number of octs
        if curdom > self.max_domain: return 0
        cdef OctAllocationContainer *cont = self.domains[curdom - 1]
        cdef int initial = cont.n_assigned
        cdef int in_boundary = 0
        # How do we bootstrap ourselves?
        for p in range(no):
            #for every oct we're trying to add find the 
            #floating point unitary position on this level
            in_boundary = 0
            for i in range(3):
                pp[i] = pos[p, i]
                dds[i] = (self.DRE[i] - self.DLE[i])/self.nn[i]
                ind[i] = <np.int64_t> ((pp[i] - self.DLE[i])/dds[i])
                cp[i] = (ind[i] + 0.5) * dds[i] + self.DLE[i]
                if ind[i] < 0 or ind[i] >= self.nn[i]:
                    in_boundary = 1
            if skip_boundary == in_boundary == 1: continue
            cur = self.next_root(curdom, ind)
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
                cur = self.next_child(curdom, ind, cur)
            # Now we should be at the right level
            cur.domain = curdom
            if local == 1:
                cur.file_ind = p
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
        ci = 0
        for oi in range(cur.n_assigned):
            o = &cur.my_octs[oi]
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        ii = ((k*2)+j)*2+i
                        if mask[o.domain_ind, ii] == 0: continue
                        coords[ci, 0] = (o.pos[0] << 1) + i
                        coords[ci, 1] = (o.pos[1] << 1) + j
                        coords[ci, 2] = (o.pos[2] << 1) + k
                        ci += 1
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
                levels[ci] = o.level
                ci += 1
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
        for oi in range(cur.n_assigned):
            o = &cur.my_octs[oi]
            for i in range(8):
                if mask[o.domain_ind, i] == 0: continue
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
        ci = 0
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
                        if mask[o.domain_ind, ii] == 0: continue
                        coords[ci, 0] = pos[0] + dx[0] * i
                        coords[ci, 1] = pos[1] + dx[1] * j
                        coords[ci, 2] = pos[2] + dx[2] * k
                        ci += 1
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
                for i in range(2):
                    for j in range(2):
                        for k in range(2):
                            ii = ((k*2)+j)*2+i
                            if mask[o.domain_ind, ii] == 0: continue
                            if o.level == level:
                                dest[local_filled] = \
                                    source[o.file_ind, ii]
                            local_filled += 1
        return local_filled

cdef class ARTOctreeContainer(RAMSESOctreeContainer):

    @cython.boundscheck(True)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_level(self, int domain, int level, dest_fields, source_fields,
                   np.ndarray[np.uint8_t, ndim=2, cast=True] mask, int offset,
                   np.int64_t subchunk_offset, np.int64_t subchunk_max):
        #Only minorly different from the RAMSES version
        #The source array is in chunks, just stop when we hit the end
        cdef np.ndarray[np.float64_t, ndim=2] source
        cdef np.ndarray[np.float64_t, ndim=1] dest
        cdef OctAllocationContainer *dom = self.domains[domain - 1]
        cdef Oct *o
        cdef int n
        cdef int i, j, k, ii
        cdef int local_pos, local_filled
        cdef np.float64_t val
        cdef np.int64_t index
        for key in dest_fields:
            local_filled = 0
            dest = dest_fields[key]
            source = source_fields[key]
            for n in range(dom.n):
                o = &dom.my_octs[n]
                index = o.file_ind-subchunk_offset
                if o.level != level: continue
                if index < 0: continue
                if index >= subchunk_max: 
                    #if we hit the end of the array,
                    #immeditely discontinue
                    return local_filled
                for i in range(2):
                    for j in range(2):
                        for k in range(2):
                            ii = ((k*2)+j)*2+i
                            if mask[o.domain_ind, ii] == 0: continue
                            dest[local_filled + offset] = \
                                source[index,ii]
                            local_filled += 1
        return local_filled


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_level_from_grid(self, int domain, int level, dest_fields, 
                             source_fields, 
                             np.ndarray[np.uint8_t, ndim=2, cast=True] mask,
                             int offset):
        #Fill  level, but instead of assuming that the source
        #order is that of the oct order, we look up the oct position
        #and fill its children from the the source field
        #As a result, source is 3D grid with 8 times as many
        #elements as the number of octs on this level in this domain
        #and with the shape of an equal-sided cube
        cdef np.ndarray[np.float64_t, ndim=3] source
        cdef np.ndarray[np.float64_t, ndim=1] dest
        cdef OctAllocationContainer *dom = self.domains[domain - 1]
        cdef Oct *o
        cdef int n
        cdef int i, j, k, ii
        cdef int local_pos, local_filled
        cdef np.float64_t val
        cdef np.int64_t ox,oy,oz
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
                            if mask[o.domain_ind, ii] == 0: continue
                            ox = (o.pos[0] << 1) + i
                            oy = (o.pos[1] << 1) + j
                            oz = (o.pos[2] << 1) + k
                            dest[local_filled + offset] = source[ox,oy,oz]
                            local_filled += 1
        return local_filled

# Now some visitor functions

cdef void visit_count_octs(Oct *o, OctVisitorData *data):
    # Number of cells visited
    data.index += 1

cdef void visit_count_total_octs(Oct *o, OctVisitorData *data):
    # Number of *octs* visited.
    if data.last != o.domain_ind:
        data.index += 1
        data.last = o.domain_ind

cdef void visit_mark_octs(Oct *o, OctVisitorData *data):
    cdef int i
    cdef np.uint8_t *arr = <np.uint8_t *> data.array
    if data.last != o.domain_ind:
        data.last = o.domain_ind
        data.index += 1
    cdef np.int64_t index = data.index * 8
    index += ((data.ind[2]*2)+data.ind[1])*2+data.ind[0] 
    arr[index] = 1

cdef void visit_index_octs(Oct *o, OctVisitorData *data):
    cdef int i
    cdef np.int64_t *arr
    if data.last != o.domain_ind:
        data.last = o.domain_ind
        arr = <np.int64_t *> data.array
        arr[o.domain_ind] = data.index
        data.index += 1

cdef void visit_icoords_octs(Oct *o, OctVisitorData *data):
    cdef np.int64_t *coords = <np.int64_t*> data.array
    cdef int i
    for i in range(3):
        coords[data.index * 3 + i] = (o.pos[i] << 1) + data.ind[i]
    data.index += 1

cdef void visit_ires_octs(Oct *o, OctVisitorData *data):
    cdef np.int64_t *ires = <np.int64_t*> data.array
    ires[data.index] = o.level
    data.index += 1

cdef void visit_fcoords_octs(Oct *o, OctVisitorData *data):
    # Note that this does not actually give the correct floating point
    # coordinates.  It gives them in some unit system where the domain is 1.0
    # in all directions, and assumes that they will be scaled later.
    cdef np.float64_t *fcoords = <np.float64_t*> data.array
    cdef int i
    cdef np.float64_t c, dx 
    dx = 1.0 / (2 << o.level)
    for i in range(3):
        c = <np.float64_t> ((o.pos[i] << 1 ) + data.ind[i]) 
        fcoords[data.index * 3 + i] = (c + 0.5) * dx
    data.index += 1

cdef void visit_fwidth_octs(Oct *o, OctVisitorData *data):
    # Note that this does not actually give the correct floating point
    # coordinates.  It gives them in some unit system where the domain is 1.0
    # in all directions, and assumes that they will be scaled later.
    cdef np.float64_t *fwidth = <np.float64_t*> data.array
    cdef int i
    cdef np.float64_t dx 
    dx = 1.0 / (2 << o.level)
    for i in range(3):
        fwidth[data.index * 3 + i] = dx
    data.index += 1
