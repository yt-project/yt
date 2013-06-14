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

from oct_container cimport OctreeContainer, Oct, OctInfo
from libc.stdlib cimport malloc, free, qsort
from libc.math cimport floor
from fp_utils cimport *
cimport numpy as np
import numpy as np
from oct_container cimport Oct, OctAllocationContainer, \
    OctreeContainer, ORDER_MAX
from selection_routines cimport SelectorObject
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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def icoords(self, int domain_id,
                np.ndarray[np.uint8_t, ndim=2, cast=True] mask,
                np.int64_t cell_count,
                np.ndarray[np.int64_t, ndim=1] level_counts):
        #Return the integer positions of the cells
        #Limited to this domain and within the mask
        #Positions are binary; aside from the root mesh
        #to each digit we just add a << 1 and a 0 or 1 
        #for each child recursively
        cdef np.ndarray[np.int64_t, ndim=2] coords
        coords = np.empty((cell_count, 3), dtype="int64")
        cdef int oi, i, ci, ii
        ci = 0
        for oi in range(self.nocts):
            o = self.oct_list[oi]
            if o.domain != domain_id: continue
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        ii = ((k*2)+j)*2+i
                        if mask[oi, ii] == 1:
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
        #Return the 'resolution' of each cell; ie the level
        cdef np.ndarray[np.int64_t, ndim=1] res
        res = np.empty(cell_count, dtype="int64")
        cdef int oi, i, ci
        ci = 0
        for oi in range(self.nocts):
            o = self.oct_list[oi]
            if o.domain != domain_id: continue
            for i in range(8):
                if mask[oi, i] == 1:
                    res[ci] = o.level
                    ci += 1
        return res

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fcoords(self, int domain_id,
                np.ndarray[np.uint8_t, ndim=2, cast=True] mask,
                np.int64_t cell_count,
                np.ndarray[np.int64_t, ndim=1] level_counts):
        #Return the floating point unitary position of every cell
        cdef np.ndarray[np.float64_t, ndim=2] coords
        coords = np.empty((cell_count, 3), dtype="float64")
        cdef int oi, i, ci
        cdef np.float64_t base_dx[3], dx[3], pos[3]
        for i in range(3):
            # This is the base_dx, but not the base distance from the center
            # position.  Note that the positions will also all be offset by
            # dx/2.0.  This is also for *oct grids*, not cells.
            base_dx[i] = (self.DRE[i] - self.DLE[i])/self.nn[i]
        ci = 0
        cdef int proc
        for oi in range(self.nocts):
            proc = 0
            for i in range(8):
                if mask[oi, i] == 1:
                    proc = 1
                    break
            if proc == 0: continue
            o = self.oct_list[oi]
            if o.domain != domain_id: continue
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
                        if mask[oi, ii] == 0: continue
                        coords[ci, 0] = pos[0] + dx[0] * i
                        coords[ci, 1] = pos[1] + dx[1] * j
                        coords[ci, 2] = pos[2] + dx[2] * k
                        ci += 1
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
        self.max_level = max_level
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
    def count_levels(self, int max_level, int domain_id,
                     np.ndarray[np.uint8_t, ndim=2, cast=True] mask):
        cdef np.ndarray[np.int64_t, ndim=1] level_count
        cdef Oct *o
        cdef int oi, i
        level_count = np.zeros(max_level+1, 'int64')
        cdef np.int64_t ndo, doff
        ndo = self.dom_offsets[domain_id + 2] \
            - self.dom_offsets[domain_id + 1]
        doff = self.dom_offsets[domain_id + 1]
        for oi in range(ndo):
            o = self.oct_list[oi + doff]
            for i in range(8):
                if mask[o.domain_ind, i] == 0: continue
                level_count[o.level] += 1
        return level_count

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

    def domain_identify(self, np.ndarray[np.uint8_t, ndim=2, cast=True] mask):
        #Return an array of length # of domains
        #Every element is True if there is at least one
        #fully refined *cell* in that domain that isn't masked out
        cdef int i, oi, m
        cdef Oct *o
        cdef np.ndarray[np.uint8_t, ndim=1, cast=True] dmask
        dmask = np.zeros(self.max_domain+1, dtype='uint8')
        for oi in range(self.nocts):
            m = 0
            o = self.oct_list[oi]
            #if o.sd.np <= 0 or o.domain == -1: continue
            for i in range(8):
                if mask[oi, i] == 1:
                    m = 1
                    break
            if m == 0: continue
            dmask[o.domain] = 1
        return dmask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_cells(self, SelectorObject selector,
              np.ndarray[np.uint8_t, ndim=2, cast=True] mask):
        #Count how many cells per level there are
        cdef int i, j, k, oi
        # pos here is CELL center, not OCT center.
        cdef np.float64_t pos[3]
        cdef int n = mask.shape[0]
        cdef int eterm[3]
        cdef np.ndarray[np.int64_t, ndim=1] count
        count = np.zeros(self.max_domain + 1, 'int64')
        for oi in range(n):
            o = self.oct_list[oi]
            if o.domain == -1: continue
            for i in range(8):
                count[o.domain] += mask[oi,i]
        return count

    def domain_and(self, np.ndarray[np.uint8_t, ndim=2, cast=True] mask,
                   int domain_id):
        cdef np.int64_t i, oi, n, use
        cdef Oct *o
        cdef np.ndarray[np.uint8_t, ndim=2] m2 = \
                np.zeros((mask.shape[0], 8), 'uint8')
        n = mask.shape[0]
        for oi in range(n):
            o = self.oct_list[oi]
            if o.domain != domain_id: continue
            use = 0
            for i in range(8):
                m2[o.domain_ind, i] = mask[o.domain_ind, i]
        return m2

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
        cdef Oct *o
        n = mask.shape[0]
        nm = 0
        # This could perhaps be faster if we 
        for oi in range(n):
            o = self.oct_list[oi]
            if o.domain != domain_id: continue
            use = 0
            for i in range(8):
                if mask[o.domain_ind, i] == 1: use = 1
            nm += use
        cdef np.ndarray[np.uint8_t, ndim=4] m2 = \
                np.zeros((2, 2, 2, nm), 'uint8')
        nm = 0
        for oi in range(n):
            o = self.oct_list[oi]
            if o.domain != domain_id: continue
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
        # Here we once again do something similar to the other functions.  We
        # need a set of indices into the final reduced, masked values.  The
        # indices will be domain.n long, and will be of type int64.  This way,
        # we can get the Oct through a .get() call, then use Oct.file_ind as an
        # index into this newly created array, then finally use the returned
        # index into the domain subset array for deposition.
        cdef np.int64_t i, j, k, oi, noct, n, nm, use, offset
        cdef Oct *o
        # For particle octrees, domain 0 is special and means non-leaf nodes.
        offset = self.dom_offsets[domain_id + 1]
        noct = self.dom_offsets[domain_id + 2] - offset
        cdef np.ndarray[np.int64_t, ndim=1] ind = np.zeros(noct, 'int64')
        nm = 0
        for oi in range(noct):
            ind[oi] = -1
            o = self.oct_list[oi + offset]
            use = 0
            for i in range(8):
                if mask[o.domain_ind, i] == 1: use = 1
            if use == 1:
                ind[oi] = nm
            nm += use
        return ind

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
