"""
Oct container tuned for Particles

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Author: Christopher Moody <chris.e.moody@gmail.com>
Affiliation: UC Santa Cruz
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2013 Matthew Turk.  All Rights Reserved.

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

from oct_container cimport OctreeContainer, Oct, OctInfo, ORDER_MAX, \
    visit_icoords_octs, visit_ires_octs, \
    visit_fcoords_octs, visit_fwidth_octs
from libc.stdlib cimport malloc, free, qsort
from libc.math cimport floor
from fp_utils cimport *
cimport numpy as np
import numpy as np
from selection_routines cimport SelectorObject, \
    OctVisitorData, oct_visitor_function
cimport cython

cdef class ParticleOctreeContainer(OctreeContainer):
    cdef Oct** oct_list
    #The starting oct index of each domain
    cdef np.int64_t *dom_offsets 
    cdef public int max_level
    #How many particles do we keep befor refining
    cdef public int n_ref

    def allocate_root(self):
        cdef int i, j, k
        cdef Oct *cur
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                for k in range(self.nn[2]):
                    cur = self.allocate_oct()
                    cur.level = 0
                    cur.pos[0] = i
                    cur.pos[1] = j
                    cur.pos[2] = k
                    cur.parent = NULL
                    self.root_mesh[i][j][k] = cur

    def __dealloc__(self):
        #Call the freemem ops on every ocy
        #of the root mesh recursively
        cdef i, j, k
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                for k in range(self.nn[2]):
                    self.visit_free(self.root_mesh[i][j][k])
        free(self.oct_list)
        free(self.dom_offsets)

    cdef void visit_free(self, Oct *o):
        #Free the memory for this oct recursively
        cdef int i, j, k
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if o.children[i][j][k] == NULL: continue
                    self.visit_free(o.children[i][j][k])
        free(o)

    def clear_fileind(self):
        cdef i, j, k
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                for k in range(self.nn[2]):
                    self.visit_clear(self.root_mesh[i][j][k])

    cdef void visit_clear(self, Oct *o):
        #Free the memory for this oct recursively
        cdef int i, j, k
        o.file_ind = 0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if o.children[i][j][k] == NULL: continue
                    self.visit_clear(o.children[i][j][k])

    def __iter__(self):
        #Get the next oct, will traverse domains
        #Note that oct containers can be sorted 
        #so that consecutive octs are on the same domain
        cdef int oi
        cdef Oct *o
        for oi in range(self.nocts):
            o = self.oct_list[oi]
            yield (o.file_ind, o.domain_ind, o.domain)

    #@cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def icoords(self, SelectorObject selector, np.uint64_t num_cells = -1):
        if num_cells == -1:
            num_cells = selector.count_octs(self)
        cdef np.ndarray[np.int64_t, ndim=2] coords
        coords = np.empty((num_cells, 3), dtype="int64")
        cdef OctVisitorData data
        data.array = <void *> coords.data
        data.index = 0
        self.visit_all_octs(selector, visit_icoords_octs, &data)
        return coords

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def ires(self, SelectorObject selector, np.uint64_t num_cells = -1):
        if num_cells == -1:
            num_cells = selector.count_octs(self)
        #Return the 'resolution' of each cell; ie the level
        cdef np.ndarray[np.int64_t, ndim=1] res
        res = np.empty(num_cells, dtype="int64")
        cdef OctVisitorData data
        data.array = <void *> res.data
        data.index = 0
        self.visit_all_octs(selector, visit_ires_octs, &data)
        return res

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fwidth(self, SelectorObject selector, np.uint64_t num_cells = -1):
        if num_cells == -1:
            num_cells = selector.count_octs(self)
        cdef np.ndarray[np.float64_t, ndim=2] fwidth
        fwidth = np.empty((num_cells, 3), dtype="float64")
        cdef OctVisitorData data
        data.array = <void *> fwidth.data
        data.index = 0
        self.visit_all_octs(selector, visit_fwidth_octs, &data)
        cdef np.float64_t base_dx
        for i in range(3):
            base_dx = (self.DRE[i] - self.DLE[i])/self.nn[i]
            fwidth[:,i] *= base_dx
        return fwidth

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fcoords(self, SelectorObject selector, np.uint64_t num_cells = -1):
        if num_cells == -1:
            num_cells = selector.count_octs(self)
        #Return the floating point unitary position of every cell
        cdef np.ndarray[np.float64_t, ndim=2] coords
        coords = np.empty((num_cells, 3), dtype="float64")
        cdef OctVisitorData data
        data.array = <void *> coords.data
        data.index = 0
        self.visit_all_octs(selector, visit_fcoords_octs, &data)
        cdef int i
        cdef np.float64_t base_dx
        for i in range(3):
            base_dx = (self.DRE[i] - self.DLE[i])/self.nn[i]
            coords[:,i] *= base_dx
            coords[:,i] += self.DLE[i]
        return coords

    def allocate_domains(self, domain_counts):
        pass

    def finalize(self):
        #This will sort the octs in the oct list
        #so that domains appear consecutively
        #And then find the oct index/offset for
        #every domain
        cdef int max_level = 0
        self.oct_list = <Oct**> malloc(sizeof(Oct*)*self.nocts)
        cdef np.int64_t i = 0, lpos = 0
        cdef int cur_dom = -1
        # We always need at least 2, and if max_domain is 0, we need 3.
        for i in range(self.nn[0]):
            for j in range(self.nn[1]):
                for k in range(self.nn[2]):
                    self.visit_assign(self.root_mesh[i][j][k], &lpos)
        assert(lpos == self.nocts)
        for i in range(self.nocts):
            self.oct_list[i].domain_ind = i
            self.oct_list[i].file_ind = -1
            max_level = imax(max_level, self.oct_list[i].level)
        self.max_level = max_level

    cdef visit_assign(self, Oct *o, np.int64_t *lpos):
        cdef int i, j, k
        self.oct_list[lpos[0]] = o
        lpos[0] += 1
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if o.children[i][j][k] != NULL:
                        self.visit_assign(o.children[i][j][k], lpos)
        return

    cdef np.int64_t get_domain_offset(self, int domain_id):
        return self.dom_offsets[domain_id + 1]

    cdef Oct* allocate_oct(self):
        #Allocate the memory, set to NULL or -1
        #We reserve space for n_ref particles, but keep
        #track of how many are used with np initially 0
        self.nocts += 1
        cdef Oct *my_oct = <Oct*> malloc(sizeof(Oct))
        cdef int i, j, k
        my_oct.domain = -1
        my_oct.file_ind = 0
        my_oct.domain_ind = self.nocts - 1
        my_oct.pos[0] = my_oct.pos[1] = my_oct.pos[2] = -1
        my_oct.level = -1
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    my_oct.children[i][j][k] = NULL
        my_oct.parent = NULL
        return my_oct

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def add(self, np.ndarray[np.uint64_t, ndim=1] indices):
        #Add this particle to the root oct
        #Then if that oct has children, add it to them recursively
        #If the child needs to be refined because of max particles, do so
        cdef np.int64_t no = indices.shape[0], p, index
        cdef int i, level, ind[3]
        if self.root_mesh[0][0][0] == NULL: self.allocate_root()
        cdef np.uint64_t *data = <np.uint64_t *> indices.data
        for p in range(no):
            # We have morton indices, which means we choose left and right by
            # looking at (MAX_ORDER - level) & with the values 1, 2, 4.
            level = 0
            index = indices[p]
            for i in range(3):
                ind[i] = (index >> ((ORDER_MAX - level)*3 + (2 - i))) & 1
            cur = self.root_mesh[ind[0]][ind[1]][ind[2]]
            if cur == NULL:
                raise RuntimeError
            while (cur.file_ind + 1) > self.n_ref:
                if level >= ORDER_MAX: break # Just dump it here.
                level += 1
                for i in range(3):
                    ind[i] = (index >> ((ORDER_MAX - level)*3 + (2 - i))) & 1
                if cur.children[ind[0]][ind[1]][ind[2]] == NULL:
                    cur = self.refine_oct(cur, index)
                    self.filter_particles(cur, data, p)
                else:
                    cur = cur.children[ind[0]][ind[1]][ind[2]]
            cur.file_ind += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Oct *refine_oct(self, Oct *o, np.uint64_t index):
        #Allocate and initialize child octs
        #Attach particles to child octs
        #Remove particles from this oct entirely
        cdef int i, j, k, m, n, ind[3]
        cdef Oct *noct
        cdef np.uint64_t prefix1, prefix2
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    noct = self.allocate_oct()
                    noct.domain = o.domain
                    noct.file_ind = 0
                    noct.level = o.level + 1
                    noct.pos[0] = (o.pos[0] << 1) + i
                    noct.pos[1] = (o.pos[1] << 1) + j
                    noct.pos[2] = (o.pos[2] << 1) + k
                    noct.parent = o
                    o.children[i][j][k] = noct
        o.file_ind = self.n_ref + 1
        for i in range(3):
            ind[i] = (index >> ((ORDER_MAX - (o.level + 1))*3 + (2 - i))) & 1
        noct = o.children[ind[0]][ind[1]][ind[2]]
        return noct

    cdef void filter_particles(self, Oct *o, np.uint64_t *data, np.int64_t p):
        # Now we look at the last nref particles to decide where they go.
        cdef int n = imin(p, self.n_ref)
        cdef np.uint64_t *arr = data + imax(p - self.n_ref, 0)
        # Now we figure out our prefix, which is the oct address at this level.
        # As long as we're actually in Morton order, we do not need to worry
        # about *any* of the other children of the oct.
        prefix1 = data[p] >> (ORDER_MAX - o.level)*3
        for i in range(n):
            prefix2 = arr[i] >> (ORDER_MAX - o.level)*3
            if (prefix1 == prefix2):
                o.file_ind += 1
        #print ind[0], ind[1], ind[2], o.file_ind, o.level

    def recursively_count(self):
        #Visit every cell, accumulate the # of cells per level
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

cdef class ParticleRegions:
    cdef np.float64_t left_edge[3]
    cdef np.float64_t dds[3]
    cdef np.float64_t idds[3]
    cdef np.int32_t dims[3]
    cdef public int nfiles
    cdef public object masks

    def __init__(self, left_edge, right_edge, dims, nfiles):
        cdef int i
        self.nfiles = nfiles
        for i in range(3):
            self.left_edge[i] = left_edge[i]
            self.dims[i] = dims[i]
            self.dds[i] = (right_edge[i] - left_edge[i])/dims[i]
            self.idds[i] = 1.0/self.dds[i]
        # We use 64-bit masks
        self.masks = []
        for i in range(nfiles/64 + 1):
            self.masks.append(np.zeros(dims, dtype="uint64"))

    def add_data_file(self, np.ndarray[np.float64_t, ndim=2] pos, int file_id):
        cdef np.int64_t no = pos.shape[0]
        cdef np.int64_t p
        cdef int ind[3], i
        cdef np.ndarray[np.uint64_t, ndim=3] mask
        mask = self.masks[file_id/64]
        val = 1 << (file_id - (file_id/64)*64)
        for p in range(no):
            # Now we locate the particle
            for i in range(3):
                ind[i] = <int> ((pos[p, i] - self.left_edge[i])*self.idds[i])
            mask[ind[0],ind[1],ind[2]] |= val

    def identify_data_files(self, SelectorObject selector):
        # This is relatively cheap to iterate over.
        cdef int i, j, k, n
        cdef np.uint64_t fmask, offset
        cdef np.float64_t LE[3], RE[3]
        cdef np.ndarray[np.uint64_t, ndim=3] mask
        files = []
        for n in range(len(self.masks)):
            fmask = 0
            mask = self.masks[n]
            LE[0] = self.left_edge[0]
            RE[0] = LE[0] + self.dds[0]
            for i in range(self.dims[0]):
                LE[1] = self.left_edge[1]
                RE[1] = LE[1] + self.dds[1]
                for j in range(self.dims[1]):
                    LE[2] = self.left_edge[2]
                    RE[2] = LE[2] + self.dds[2]
                    for k in range(self.dims[2]):
                        if selector.select_grid(LE, RE, 0) == 1:
                            fmask |= mask[i,j,k]
                        LE[2] += self.dds[2]
                        RE[2] += self.dds[2]
                    LE[1] += self.dds[1]
                    RE[1] += self.dds[1]
                LE[0] += self.dds[0]
                RE[0] += self.dds[0]
            # Now we iterate through...
            for i in range(64):
                if ((fmask >> i) & 1) == 1:
                    files.append(i + n * 64)
        return files
