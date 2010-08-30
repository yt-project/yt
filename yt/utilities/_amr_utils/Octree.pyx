"""
A refine-by-two AMR-specific octree

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

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


import numpy as np
cimport numpy as np
# Double up here for def'd functions
cimport numpy as cnp
cimport cython

from stdlib cimport malloc, free, abs

cdef extern from "stdlib.h":
    # NOTE that size_t might not be int
    void *alloca(int)

cdef struct OctreeNode:
    np.float64_t *val
    np.float64_t weight_val
    np.int64_t pos[3]
    int level
    int nvals
    OctreeNode *children[2][2][2]

cdef void OTN_add_value(OctreeNode *self,
        np.float64_t *val, np.float64_t weight_val):
    cdef int i
    for i in range(self.nvals):
        self.val[i] += val[i]
    self.weight_val += weight_val

cdef void OTN_refine(OctreeNode *self):
    cdef int i, j, i1, j1
    cdef np.int64_t npos[3]
    cdef OctreeNode *node
    for i in range(2):
        npos[0] = self.pos[0] * 2 + i
        for j in range(2):
            npos[1] = self.pos[1] * 2 + j
            # We have to be careful with allocation...
            for k in range(2):
                npos[2] = self.pos[2] * 2 + k
                self.children[i][j][k] = OTN_initialize(
                            npos,
                            self.nvals, self.val, self.weight_val,
                            self.level + 1)
    for i in range(self.nvals): self.val[i] = 0.0
    self.weight_val = 0.0

cdef OctreeNode *OTN_initialize(np.int64_t pos[3], int nvals,
                        np.float64_t *val, np.float64_t weight_val,
                        int level):
    cdef OctreeNode *node
    cdef int i, j
    node = <OctreeNode *> malloc(sizeof(OctreeNode))
    node.pos[0] = pos[0]
    node.pos[1] = pos[1]
    node.pos[2] = pos[2]
    node.nvals = nvals
    node.val = <np.float64_t *> malloc(
                nvals * sizeof(np.float64_t))
    for i in range(nvals):
        node.val[i] = val[i]
    node.weight_val = weight_val
    for i in range(2):
        for j in range(2):
            for k in range(2):
                node.children[i][j][k] = NULL
    node.level = level
    return node

cdef void OTN_free(OctreeNode *node):
    cdef int i, j
    for i in range(2):
        for j in range(2):
            for k in range(2):
                if node.children[i][j][k] == NULL: continue
                OTN_free(node.children[i][j][k])
    free(node.val)
    free(node)

cdef class Octree:
    cdef int nvals
    cdef np.int64_t po2[80]
    cdef OctreeNode ****root_nodes
    cdef np.int64_t top_grid_dims[3]

    def __cinit__(self, np.ndarray[np.int64_t, ndim=1] top_grid_dims,
                  int nvals):
        cdef int i, j
        cdef OctreeNode *node
        cdef np.int64_t pos[3]
        cdef np.float64_t *vals = <np.float64_t *> alloca(
                sizeof(np.float64_t)*nvals)
        cdef np.float64_t weight_val = 0.0
        self.nvals = nvals
        for i in range(nvals): vals[i] = 0.0

        self.top_grid_dims[0] = top_grid_dims[0]
        self.top_grid_dims[1] = top_grid_dims[1]
        self.top_grid_dims[2] = top_grid_dims[2]

        # This wouldn't be necessary if we did bitshifting...
        for i in range(80):
            self.po2[i] = 2**i
        # Cython doesn't seem to like sizeof(OctreeNode ***)
        self.root_nodes = <OctreeNode ****> \
            malloc(sizeof(void*) * top_grid_dims[0])

        # We initialize our root values to 0.0.
        for i in range(top_grid_dims[0]):
            pos[0] = i
            self.root_nodes[i] = <OctreeNode ***> \
                malloc(sizeof(OctreeNode **) * top_grid_dims[1])
            for j in range(top_grid_dims[1]):
                pos[1] = j
                self.root_nodes[i][j] = <OctreeNode **> \
                    malloc(sizeof(OctreeNode *) * top_grid_dims[1])
                for k in range(top_grid_dims[2]):
                    pos[2] = k
                    self.root_nodes[i][j][k] = OTN_initialize(
                        pos, nvals, vals, weight_val, 0)

    cdef void add_to_position(self,
                 int level, np.int64_t pos[3],
                 np.float64_t *val,
                 np.float64_t weight_val):
        cdef int i, j
        cdef OctreeNode *node
        node = self.find_on_root_level(pos, level)
        cdef np.int64_t fac
        for L in range(level):
            if node.children[0][0][0] == NULL:
                OTN_refine(node)
            # Maybe we should use bitwise operators?
            fac = self.po2[level - L - 1]
            i = (pos[0] >= fac*(2*node.pos[0]+1))
            j = (pos[1] >= fac*(2*node.pos[1]+1))
            k = (pos[2] >= fac*(2*node.pos[2]+1))
            node = node.children[i][j][k]
        OTN_add_value(node, val, weight_val)
            
    cdef OctreeNode *find_on_root_level(self, np.int64_t pos[3], int level):
        # We need this because the root level won't just have four children
        # So we find on the root level, then we traverse the tree.
        cdef np.int64_t i, j
        i = <np.int64_t> (pos[0] / self.po2[level])
        j = <np.int64_t> (pos[1] / self.po2[level])
        k = <np.int64_t> (pos[2] / self.po2[level])
        return self.root_nodes[i][j][k]
        
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def add_array_to_tree(self, int level,
            np.ndarray[np.int64_t, ndim=1] pxs,
            np.ndarray[np.int64_t, ndim=1] pys,
            np.ndarray[np.int64_t, ndim=1] pzs,
            np.ndarray[np.float64_t, ndim=2] pvals,
            np.ndarray[np.float64_t, ndim=1] pweight_vals):
        cdef int np = pxs.shape[0]
        cdef int p
        cdef cnp.float64_t *vals
        cdef cnp.float64_t *data = <cnp.float64_t *> pvals.data
        cdef cnp.int64_t pos[3]
        for p in range(np):
            vals = data + self.nvals*p
            pos[0] = pxs[p]
            pos[1] = pys[p]
            pos[2] = pzs[p]
            self.add_to_position(level, pos, vals, pweight_vals[p])

    def add_grid_to_tree(self, int level,
                         np.ndarray[np.int64_t, ndim=1] start_index,
                         np.ndarray[np.float64_t, ndim=2] pvals,
                         np.ndarray[np.float64_t, ndim=2] wvals,
                         np.ndarray[np.int32_t, ndim=2] cm):
        pass

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def get_all_from_level(self, int level, int count_only = 0):
        cdef int i, j
        cdef int total = 0
        vals = []
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                for k in range(self.top_grid_dims[2]):
                    total += self.count_at_level(self.root_nodes[i][j][k], level)
        if count_only: return total
        # Allocate our array
        cdef np.ndarray[np.int64_t, ndim=2] npos
        cdef np.ndarray[np.float64_t, ndim=2] nvals
        cdef np.ndarray[np.float64_t, ndim=1] nwvals
        npos = np.zeros( (total, 2), dtype='int64')
        nvals = np.zeros( (total, self.nvals), dtype='float64')
        nwvals = np.zeros( total, dtype='float64')
        cdef np.int64_t curpos = 0
        cdef np.int64_t *pdata = <np.int64_t *> npos.data
        cdef np.float64_t *vdata = <np.float64_t *> nvals.data
        cdef np.float64_t *wdata = <np.float64_t *> nwvals.data
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                for k in range(self.top_grid_dims[2]):
                    curpos += self.fill_from_level(self.root_nodes[i][j][k],
                        level, curpos, pdata, vdata, wdata)
        return npos, nvals, nwvals

    cdef int count_at_level(self, OctreeNode *node, int level):
        cdef int i, j
        # We only really return a non-zero, calculated value if we are at the
        # level in question.
        if node.level == level:
            # We return 1 if there are no finer points at this level and zero
            # if there are
            return (node.children[0][0][0] == NULL)
        if node.children[0][0][0] == NULL: return 0
        cdef int count = 0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    count += self.count_at_level(node.children[i][j][k], level)
        return count

    cdef int fill_from_level(self, OctreeNode *node, int level,
                              np.int64_t curpos,
                              np.int64_t *pdata,
                              np.float64_t *vdata,
                              np.float64_t *wdata):
        cdef int i, j
        if node.level == level:
            if node.children[0][0][0] != NULL: return 0
            for i in range(self.nvals):
                vdata[self.nvals * curpos + i] = node.val[i]
            wdata[curpos] = node.weight_val
            pdata[curpos * 3] = node.pos[0]
            pdata[curpos * 3 + 1] = node.pos[1]
            pdata[curpos * 3 + 2] = node.pos[2]
            return 1
        if node.children[0][0] == NULL: return 0
        cdef np.int64_t added = 0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    added += self.fill_from_level(node.children[i][j][k],
                            level, curpos + added, pdata, vdata, wdata)
        return added

    def __dealloc__(self):
        cdef int i, j
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                for k in range(self.top_grid_dims[2]):
                    OTN_free(self.root_nodes[i][j][k])
                free(self.root_nodes[i][j])
            free(self.root_nodes[i])
        free(self.root_nodes)

