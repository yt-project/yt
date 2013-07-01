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
import selection_routines
cimport oct_visitors
from oct_visitors cimport cind
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
        oct.file_ind = oct.domain = -1
        oct.domain_ind = n + n_cont.offset
        oct.children = NULL
    if prev != NULL:
        prev.next = n_cont
    n_cont.next = NULL
    return n_cont

cdef void free_octs(
        OctAllocationContainer *first):
    cdef OctAllocationContainer *cur
    while first != NULL:
        cur = first
        for i in range(cur.n):
            if cur.my_octs[i].children != NULL:
                free(cur.my_octs[i].children)
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

    def __init__(self, oct_domain_dimensions, domain_left_edge,
                 domain_right_edge, partial_coverage = 0):
        # This will just initialize the root mesh octs
        self.partial_coverage = partial_coverage
        cdef int i, j, k, p
        for i in range(3):
            self.nn[i] = oct_domain_dimensions[i]
        self.max_domain = -1
        p = 0
        self.nocts = 0 # Increment when initialized
        for i in range(3):
            self.DLE[i] = domain_left_edge[i] #0
            self.DRE[i] = domain_right_edge[i] #num_grid
        self._initialize_root_mesh()

    def _initialize_root_mesh(self):
        self.root_mesh = <Oct****> malloc(sizeof(void*) * self.nn[0])
        for i in range(self.nn[0]):
            self.root_mesh[i] = <Oct ***> malloc(sizeof(void*) * self.nn[1])
            for j in range(self.nn[1]):
                self.root_mesh[i][j] = <Oct **> malloc(sizeof(void*) * self.nn[2])
                for k in range(self.nn[2]):
                    self.root_mesh[i][j][k] = NULL

    def __dealloc__(self):
        free_octs(self.cont)
        if self.root_mesh == NULL: return
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                if self.root_mesh[i][j] == NULL: continue
                free(self.root_mesh[i][j])
            if self.root_mesh[i] == NULL: continue
            free(self.root_mesh[i])
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

    @cython.cdivision(True)
    cdef void visit_all_octs(self, SelectorObject selector,
                        oct_visitor_function *func,
                        OctVisitorData *data):
        cdef int i, j, k, n, vc
        vc = self.partial_coverage
        data.global_index = -1
        data.level = 0
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
                    if self.root_mesh[i][j][k] == NULL:
                        raise RuntimeError
                    data.pos[0] = i
                    data.pos[1] = j
                    data.pos[2] = k
                    selector.recursively_visit_octs(
                        self.root_mesh[i][j][k],
                        pos, dds, 0, func, data, vc)
                    pos[2] += dds[2]
                pos[1] += dds[1]
            pos[0] += dds[0]

    cdef void oct_bounds(self, Oct *o, np.float64_t *corner, np.float64_t *size):
        cdef int i
        #for i in range(3):
        #    size[i] = (self.DRE[i] - self.DLE[i]) / (self.nn[i] << o.level)
        #    corner[i] = o.pos[i] * size[i] + self.DLE[i]

    cdef np.int64_t get_domain_offset(self, int domain_id):
        return 0

    cdef int get_root(self, int ind[3], Oct **o):
        o[0] = self.root_mesh[ind[0]][ind[1]][ind[2]]
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Oct *get(self, np.float64_t ppos[3], OctInfo *oinfo = NULL):
        #Given a floating point position, retrieve the most
        #refined oct at that time
        cdef int ind[3]
        cdef np.float64_t dds[3], cp[3], pp[3]
        cdef Oct *cur, *next
        cur = next = NULL
        cdef int i
        for i in range(3):
            dds[i] = (self.DRE[i] - self.DLE[i])/self.nn[i]
            ind[i] = <np.int64_t> ((ppos[i] - self.DLE[i])/dds[i])
            cp[i] = (ind[i] + 0.5) * dds[i] + self.DLE[i]
        self.get_root(ind, &next)
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
            if cur.children != NULL:
                next = cur.children[cind(ind[0],ind[1],ind[2])]
            else:
                next = NULL
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

    def domain_identify(self, SelectorObject selector):
        cdef np.ndarray[np.uint8_t, ndim=1] domain_mask
        domain_mask = np.zeros(self.max_domain, dtype="uint8")
        cdef OctVisitorData data
        data.array = domain_mask.data
        data.domain = -1
        self.visit_all_octs(selector, oct_visitors.identify_octs, &data)
        cdef int i
        domain_ids = []
        for i in range(self.max_domain):
            if domain_mask[i] == 1:
                domain_ids.append(i+1)
        return domain_ids

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
        raise RuntimeError
        #for ni in range(3):
        #    for nj in range(3):
        #        for nk in range(3):
        #            if ni == nj == nk == 1:
        #                neighbors[nn] = o
        #                nn += 1
        #                continue
        #            npos[0] = o.pos[0] + (ni - 1)
        #            npos[1] = o.pos[1] + (nj - 1)
        #            npos[2] = o.pos[2] + (nk - 1)
        #            for i in range(3):
        #                # Periodicity
        #                if npos[i] == -1:
        #                    npos[i] = (self.nn[i]  << o.level) - 1
        #                elif npos[i] == (self.nn[i] << o.level):
        #                    npos[i] = 0
        #                curopos[i] = o.pos[i]
        #                curnpos[i] = npos[i] 
        #            # Now we have our neighbor position and a safe place to
        #            # keep it.  curnpos will be the root index of the neighbor
        #            # at a given level, and npos will be constant.  curopos is
        #            # the candidate root at a level.
        #            candidate = o
        #            while candidate != NULL:
        #                if ((curopos[0] == curnpos[0]) and 
        #                    (curopos[1] == curnpos[1]) and
        #                    (curopos[2] == curnpos[2])):
        #                    break
        #                # This one doesn't meet it, so we pop up a level.
        #                # First we update our positions, then we update our
        #                # candidate.
        #                for i in range(3):
        #                    # We strip a digit off the right
        #                    curopos[i] = (curopos[i] >> 1)
        #                    curnpos[i] = (curnpos[i] >> 1)
        #                # Now we update to the candidate's parent, which should
        #                # have a matching position to curopos[]
        #                # TODO: This has not survived the transition to
        #                # mostly-stateless Octs!
        #                raise RuntimeError
        #                candidate = candidate.parent
        #            if candidate == NULL:
        #                # Worst case scenario
        #                for i in range(3):
        #                    ind[i] = (npos[i] >> (o.level))
        #                candidate = self.root_mesh[ind[0]][ind[1]][ind[2]]
        #            # Now we have the common root, which may be NULL
        #            while candidate.level < o.level:
        #                dl = o.level - (candidate.level + 1)
        #                for i in range(3):
        #                    ind[i] = (npos[i] >> dl) & 1
        #                if candidate.children[cind(ind[0],ind[1],ind[2])] \
        #                        == NULL:
        #                    break
        #                candidate = candidate.children[cind(ind[0],ind[1],ind[2])]
        #            neighbors[nn] = candidate
        #            nn += 1

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
        raise RuntimeError
        for i in range(27):
            self.oct_bounds(neighbors[i], corner, size)
            for ii in range(3):
                bounds[i, ii] = corner[ii]
                bounds[i, 3+ii] = size[ii]
        return bounds

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def icoords(self, SelectorObject selector, np.int64_t num_octs = -1,
                int domain_id = -1):
        if num_octs == -1:
            num_octs = selector.count_octs(self, domain_id)
        cdef np.ndarray[np.int64_t, ndim=2] coords
        coords = np.empty((num_octs * 8, 3), dtype="int64")
        cdef OctVisitorData data
        data.array = <void *> coords.data
        data.index = 0
        data.domain = domain_id
        self.visit_all_octs(selector, oct_visitors.icoords_octs, &data)
        return coords

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def ires(self, SelectorObject selector, np.int64_t num_octs = -1,
                int domain_id = -1):
        if num_octs == -1:
            num_octs = selector.count_octs(self, domain_id)
        #Return the 'resolution' of each cell; ie the level
        cdef np.ndarray[np.int64_t, ndim=1] res
        res = np.empty(num_octs * 8, dtype="int64")
        cdef OctVisitorData data
        data.array = <void *> res.data
        data.index = 0
        data.domain = domain_id
        self.visit_all_octs(selector, oct_visitors.ires_octs, &data)
        return res

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fwidth(self, SelectorObject selector, np.int64_t num_octs = -1,
                int domain_id = -1):
        if num_octs == -1:
            num_octs = selector.count_octs(self, domain_id)
        cdef np.ndarray[np.float64_t, ndim=2] fwidth
        fwidth = np.empty((num_octs * 8, 3), dtype="float64")
        cdef OctVisitorData data
        data.array = <void *> fwidth.data
        data.index = 0
        data.domain = domain_id
        self.visit_all_octs(selector, oct_visitors.fwidth_octs, &data)
        cdef np.float64_t base_dx
        for i in range(3):
            base_dx = (self.DRE[i] - self.DLE[i])/self.nn[i]
            fwidth[:,i] *= base_dx
        return fwidth

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fcoords(self, SelectorObject selector, np.int64_t num_octs = -1,
                int domain_id = -1):
        if num_octs == -1:
            num_octs = selector.count_octs(self, domain_id)
        #Return the floating point unitary position of every cell
        cdef np.ndarray[np.float64_t, ndim=2] coords
        coords = np.empty((num_octs * 8, 3), dtype="float64")
        cdef OctVisitorData data
        data.array = <void *> coords.data
        data.index = 0
        data.domain = domain_id
        self.visit_all_octs(selector, oct_visitors.fcoords_octs, &data)
        cdef int i
        cdef np.float64_t base_dx
        for i in range(3):
            base_dx = (self.DRE[i] - self.DLE[i])/self.nn[i]
            coords[:,i] *= base_dx
            coords[:,i] += self.DLE[i]
        return coords

    def selector_fill(self, SelectorObject selector,
                      np.ndarray source,
                      np.ndarray dest = None,
                      np.int64_t offset = 0, int dims = 1,
                      int domain_id = -1):
        # This is actually not correct.  The hard part is that we need to
        # iterate the same way visit_all_octs does, but we need to track the
        # number of octs total visited.
        cdef np.int64_t num_cells = -1
        if dest is None:
            # Note that RAMSES can have partial refinement inside an Oct.  This
            # means we actually do want the number of Octs, not the number of
            # cells.
            num_cells = selector.count_oct_cells(self, domain_id)
            if dims > 1:
                dest = np.zeros((num_cells, dims), dtype=source.dtype,
                    order='C')
            else:
                dest = np.zeros(num_cells, dtype=source.dtype, order='C')
        cdef OctVisitorData data
        data.index = offset
        data.domain = domain_id
        # We only need this so we can continue calculating the offset
        data.dims = dims
        cdef void *p[2]
        p[0] = source.data
        p[1] = dest.data
        data.array = &p
        cdef oct_visitor_function *func
        if source.dtype != dest.dtype:
            raise RuntimeError
        if source.dtype == np.int64:
            func = oct_visitors.copy_array_i64
        elif source.dtype == np.float64:
            func = oct_visitors.copy_array_f64
        else:
            raise NotImplementedError
        self.visit_all_octs(selector, func, &data)
        if (data.global_index + 1) * 8 * data.dims > source.size:
            print "GLOBAL INDEX RAN AHEAD.",
            print (data.global_index + 1) * 8 * data.dims - source.size
            print dest.size, source.size, num_cells
            raise RuntimeError
        if data.index > dest.size:
            print "DEST INDEX RAN AHEAD.",
            print data.index - dest.size
            raise RuntimeError
        if num_cells >= 0:
            return dest
        return data.index - offset

    def domain_ind(self, selector, int domain_id = -1):
        cdef np.ndarray[np.int64_t, ndim=1] ind
        # Here's where we grab the masked items.
        ind = np.zeros(self.nocts, 'int64') - 1
        cdef OctVisitorData data
        data.domain = domain_id
        data.array = ind.data
        data.index = 0
        data.last = -1
        self.visit_all_octs(selector, oct_visitors.index_octs, &data)
        return ind

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def add(self, int curdom, int curlevel,
            np.ndarray[np.float64_t, ndim=2] pos,
            int skip_boundary = 1):
        cdef int level, no, p, i, j, k, ind[3]
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
            if cur == NULL: raise RuntimeError
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
            cur.file_ind = p
        return cont.n_assigned - initial

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

    cdef Oct* next_root(self, int domain_id, int ind[3]):
        cdef Oct *next = self.root_mesh[ind[0]][ind[1]][ind[2]]
        if next != NULL: return next
        cdef OctAllocationContainer *cont = self.domains[domain_id - 1]
        if cont.n_assigned >= cont.n: raise RuntimeError
        next = &cont.my_octs[cont.n_assigned]
        cont.n_assigned += 1
        self.root_mesh[ind[0]][ind[1]][ind[2]] = next
        self.nocts += 1
        return next

    cdef Oct* next_child(self, int domain_id, int ind[3], Oct *parent):
        cdef int i
        cdef Oct *next = NULL
        if parent.children != NULL:
            next = parent.children[cind(ind[0],ind[1],ind[2])]
        else:
            parent.children = <Oct **> malloc(sizeof(Oct *) * 8)
            for i in range(8):
                parent.children[i] = NULL
        if next != NULL: return next
        cdef OctAllocationContainer *cont = self.domains[domain_id - 1]
        if cont.n_assigned >= cont.n: raise RuntimeError
        next = &cont.my_octs[cont.n_assigned]
        cont.n_assigned += 1
        parent.children[cind(ind[0],ind[1],ind[2])] = next
        self.nocts += 1
        return next

    def file_index_octs(self, SelectorObject selector, int domain_id,
                        num_cells = -1):
        # We create oct arrays of the correct size
        cdef np.int64_t i
        cdef np.ndarray[np.uint8_t, ndim=1] levels
        cdef np.ndarray[np.uint8_t, ndim=1] cell_inds
        cdef np.ndarray[np.int64_t, ndim=1] file_inds
        if num_cells < 0:
            num_cells = selector.count_oct_cells(self, domain_id)
        levels = np.zeros(num_cells, dtype="uint8")
        file_inds = np.zeros(num_cells, dtype="int64")
        cell_inds = np.zeros(num_cells, dtype="uint8")
        for i in range(num_cells):
            levels[i] = 100
            file_inds[i] = -1
            cell_inds[i] = 9
        cdef OctVisitorData data
        data.index = 0
        cdef void *p[3]
        p[0] = levels.data
        p[1] = file_inds.data
        p[2] = cell_inds.data
        data.array = p
        data.domain = domain_id
        self.visit_all_octs(selector, oct_visitors.fill_file_indices, &data)
        return levels, cell_inds, file_inds

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_level(self, int level,
                   np.ndarray[np.uint8_t, ndim=1] levels,
                   np.ndarray[np.uint8_t, ndim=1] cell_inds,
                   np.ndarray[np.int64_t, ndim=1] file_inds,
                   dest_fields, source_fields):
        cdef np.ndarray[np.float64_t, ndim=2] source
        cdef np.ndarray[np.float64_t, ndim=1] dest
        cdef int n
        cdef int i, di
        cdef int local_pos, local_filled
        cdef np.float64_t val
        for key in dest_fields:
            dest = dest_fields[key]
            source = source_fields[key]
            for i in range(levels.shape[0]):
                if levels[i] != level: continue
                dest[i] = source[file_inds[i], cell_inds[i]]

    def finalize(self):
        cdef SelectorObject selector = selection_routines.AlwaysSelector(None)
        cdef OctVisitorData data
        data.index = 0
        data.domain = 1
        self.visit_all_octs(selector, oct_visitors.assign_domain_ind, &data)
        assert ((data.global_index+1)*8 == data.index)

cdef int root_node_compare(void *a, void *b) nogil:
    cdef OctKey *ao, *bo
    ao = <OctKey *>a
    bo = <OctKey *>b
    if ao.key < bo.key:
        return -1
    elif ao.key == bo.key:
        return 0
    else:
        return 1

cdef class RAMSESOctreeContainer(OctreeContainer):

    def __init__(self, domain_dimensions, domain_left_edge, domain_right_edge):
        cdef int i, j, k, p
        self.partial_coverage = 1
        for i in range(3):
            self.nn[i] = domain_dimensions[i]
        self.max_domain = -1
        self.nocts = 0 # Increment when initialized
        self.root_mesh = NULL
        self.root_nodes = NULL
        self.tree_root = NULL
        self.num_root = 0
        # We don't initialize the octs yet
        for i in range(3):
            self.DLE[i] = domain_left_edge[i] #0
            self.DRE[i] = domain_right_edge[i] #num_grid

    cdef int get_root(self, int ind[3], Oct **o):
        o[0] = NULL
        cdef int i
        cdef np.int64_t key = 0
        for i in range(3):
            key |= ((<np.int64_t>ind[i]) << 20 * (2 - i))
        cdef OctKey okey, **oresult
        okey.key = key
        okey.node = NULL
        oresult = <OctKey **> tfind(<void*>&okey,
            &self.tree_root, root_node_compare)
        if oresult != NULL:
            o[0] = oresult[0].node

    @cython.cdivision(True)
    cdef void visit_all_octs(self, SelectorObject selector,
                        oct_visitor_function *func,
                        OctVisitorData *data):
        cdef int i, j, k, n, vc
        cdef np.int64_t key, ukey
        data.global_index = -1
        data.level = 0
        vc = self.partial_coverage
        cdef np.float64_t pos[3], dds[3]
        # This dds is the oct-width
        for i in range(3):
            dds[i] = (self.DRE[i] - self.DLE[i]) / self.nn[i]
        # Pos is the center of the octs
        cdef Oct *o
        ukey = 0
        for i in range(20):
            ukey |= (1 << i)
        for i in range(self.num_root):
            o = self.root_nodes[i].node
            key = self.root_nodes[i].key
            for j in range(3):
                data.pos[2 - j] = (key & ukey)
                key = key >> 20
            for j in range(3):
                pos[j] = self.DLE[j] + (data.pos[j] + 0.5) * dds[j]
            selector.recursively_visit_octs(
                o, pos, dds, 0, func, data, vc)

    cdef np.int64_t get_domain_offset(self, int domain_id):
        return 0 # We no longer have a domain offset.

    cdef Oct* next_root(self, int domain_id, int ind[3]):
        # We assume that 20 bits is enough for each index.
        cdef int i
        cdef Oct *next
        self.get_root(ind, &next)
        if next != NULL: return next
        cdef OctAllocationContainer *cont = self.domains[domain_id - 1]
        if cont.n_assigned >= cont.n:
            print "Too many assigned."
            return NULL
        if self.num_root >= self.max_root:
            print "Too many roots."
            return NULL
        next = &cont.my_octs[cont.n_assigned]
        cont.n_assigned += 1
        cdef np.int64_t key = 0
        cdef OctKey *ikey = &self.root_nodes[self.num_root]
        for i in range(3):
            key |= ((<np.int64_t>ind[i]) << 20 * (2 - i))
        self.root_nodes[self.num_root].key = key
        self.root_nodes[self.num_root].node = next
        tsearch(<void*>ikey, &self.tree_root, root_node_compare)
        self.num_root += 1
        self.nocts += 1
        return next

    def allocate_domains(self, domain_counts, int root_nodes):
        OctreeContainer.allocate_domains(self, domain_counts)
        self.root_nodes = <OctKey*> malloc(sizeof(OctKey) * root_nodes)
        self.max_root = root_nodes
        for i in range(root_nodes):
            self.root_nodes[i].key = -1
            self.root_nodes[i].node = NULL

    def __dealloc__(self):
        # This gets called BEFORE the superclass deallocation.  But, both get
        # called.
        if self.root_nodes != NULL: free(self.root_nodes)
        if self.domains != NULL: free(self.domains)

cdef class ARTOctreeContainer(OctreeContainer):

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
                #if o.level != level: continue
                for i in range(2):
                    for j in range(2):
                        for k in range(2):
                            ii = ((k*2)+j)*2+i
                            if mask[o.domain_ind, ii] == 0: continue
                            #ox = (o.pos[0] << 1) + i
                            #oy = (o.pos[1] << 1) + j
                            #oz = (o.pos[2] << 1) + k
                            dest[local_filled + offset] = source[ox,oy,oz]
                            local_filled += 1
        return local_filled
