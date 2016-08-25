"""
A refine-by-two AMR-specific quadtree



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import numpy as np
cimport numpy as np
cimport cython

from libc.stdlib cimport malloc, free, abs
from cython.operator cimport dereference as deref, preincrement as inc
from yt.utilities.lib.fp_utils cimport fmax

from yt.utilities.exceptions import YTIntDomainOverflow

cdef extern from "platform_dep.h":
    # NOTE that size_t might not be int
    void *alloca(int)

cdef struct QuadTreeNode:
    np.float64_t *val
    np.float64_t weight_val
    np.int64_t pos[2]
    QuadTreeNode *children[2][2]

ctypedef void QTN_combine(QuadTreeNode *self,
        np.float64_t *val, np.float64_t weight_val,
        int nvals)

cdef void QTN_add_value(QuadTreeNode *self,
        np.float64_t *val, np.float64_t weight_val,
        int nvals):
    cdef int i
    for i in range(nvals):
        self.val[i] += val[i]
    self.weight_val += weight_val

cdef void QTN_max_value(QuadTreeNode *self,
        np.float64_t *val, np.float64_t weight_val,
        int nvals):
    cdef int i
    for i in range(nvals):
        self.val[i] = fmax(val[i], self.val[i])
    self.weight_val = 1.0

cdef void QTN_refine(QuadTreeNode *self, int nvals):
    cdef int i, j
    cdef np.int64_t npos[2]
    cdef np.float64_t *tvals = <np.float64_t *> alloca(
            sizeof(np.float64_t) * nvals)
    for i in range(nvals): tvals[i] = 0.0
    for i in range(2):
        npos[0] = self.pos[0] * 2 + i
        for j in range(2):
            npos[1] = self.pos[1] * 2 + j
            # We have to be careful with allocation...
            self.children[i][j] = QTN_initialize(
                        npos, nvals, tvals, 0.0)

cdef QuadTreeNode *QTN_initialize(np.int64_t pos[2], int nvals,
                        np.float64_t *val, np.float64_t weight_val):
    cdef QuadTreeNode *node
    cdef int i, j
    node = <QuadTreeNode *> malloc(sizeof(QuadTreeNode))
    node.pos[0] = pos[0]
    node.pos[1] = pos[1]
    node.val = <np.float64_t *> malloc(
                nvals * sizeof(np.float64_t))
    for i in range(2):
        for j in range(2):
            node.children[i][j] = NULL
    if val != NULL:
        for i in range(nvals):
            node.val[i] = val[i]
        node.weight_val = weight_val
    return node

cdef void QTN_free(QuadTreeNode *node):
    cdef int i, j
    for i in range(2):
        for j in range(2):
            if node.children[i][j] == NULL: continue
            QTN_free(node.children[i][j])
    free(node.val)
    free(node)

cdef class QuadTree:
    cdef int nvals
    cdef QuadTreeNode ***root_nodes
    cdef np.int64_t top_grid_dims[2]
    cdef int merged
    cdef int num_cells
    cdef QTN_combine *combine
    cdef np.float64_t bounds[4]
    cdef np.float64_t dds[2]
    cdef np.int64_t last_dims[2]
    cdef int max_level

    def __cinit__(self, np.ndarray[np.int64_t, ndim=1] top_grid_dims,
                  int nvals, bounds, method = "integrate"):
        if method == "integrate":
            self.combine = QTN_add_value
        elif method == "mip":
            self.combine = QTN_max_value
        else:
            raise NotImplementedError
        self.merged = 1
        self.max_level = 0
        cdef int i, j
        cdef np.int64_t pos[2]
        cdef np.float64_t *vals = <np.float64_t *> malloc(
                sizeof(np.float64_t)*nvals)
        cdef np.float64_t weight_val = 0.0
        self.nvals = nvals
        for i in range(nvals): vals[i] = 0.0
        for i in range(4):
            self.bounds[i] = bounds[i]

        self.top_grid_dims[0] = top_grid_dims[0]
        self.top_grid_dims[1] = top_grid_dims[1]
        self.dds[0] = (self.bounds[1] - self.bounds[0])/self.top_grid_dims[0]
        self.dds[1] = (self.bounds[3] - self.bounds[2])/self.top_grid_dims[1]

        self.root_nodes = <QuadTreeNode ***> \
            malloc(sizeof(QuadTreeNode **) * top_grid_dims[0])

        # We initialize our root values to 0.0.
        for i in range(top_grid_dims[0]):
            pos[0] = i
            self.root_nodes[i] = <QuadTreeNode **> \
                malloc(sizeof(QuadTreeNode *) * top_grid_dims[1])
            for j in range(top_grid_dims[1]):
                pos[1] = j
                self.root_nodes[i][j] = QTN_initialize(
                    pos, nvals, vals, weight_val)
        self.num_cells = self.top_grid_dims[0] * self.top_grid_dims[1]
        free(vals)

    cdef int count_total_cells(self, QuadTreeNode *root):
        cdef int total = 0
        cdef int i, j
        if root.children[0][0] == NULL: return 1
        for i in range(2):
            for j in range(2):
                total += self.count_total_cells(root.children[i][j])
        return total + 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int fill_buffer(self, QuadTreeNode *root, int curpos,
                          np.ndarray[np.int32_t, ndim=1] refined,
                          np.ndarray[np.float64_t, ndim=2] values,
                          np.ndarray[np.float64_t, ndim=1] wval):
        cdef int i, j
        for i in range(self.nvals):
            values[curpos, i] = root.val[i]
        wval[curpos] = root.weight_val
        if root.children[0][0] != NULL: refined[curpos] = 1
        else: return curpos+1
        curpos += 1
        for i in range(2):
            for j in range(2):
                curpos = self.fill_buffer(root.children[i][j], curpos,
                                 refined, values, wval)
        return curpos

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int unfill_buffer(self, QuadTreeNode *root, int curpos,
                          np.ndarray[np.int32_t, ndim=1] refined,
                          np.ndarray[np.float64_t, ndim=2] values,
                          np.ndarray[np.float64_t, ndim=1] wval):
        cdef int i, j
        for i in range(self.nvals):
            root.val[i] = values[curpos, i]
        root.weight_val = wval[curpos]
        if refined[curpos] == 0: return curpos+1
        curpos += 1
        cdef QuadTreeNode *child
        cdef np.int64_t pos[2]
        for i in range(2):
            for j in range(2):
                pos[0] = root.pos[0]*2 + i
                pos[1] = root.pos[1]*2 + j
                child = QTN_initialize(pos, self.nvals, NULL, 0.0)
                root.children[i][j] = child
                curpos = self.unfill_buffer(child, curpos, refined, values, wval)
        return curpos


    @cython.boundscheck(False)
    @cython.wraparound(False)
    def frombuffer(self, np.ndarray[np.int32_t, ndim=1] refined,
                         np.ndarray[np.float64_t, ndim=2] values,
                         np.ndarray[np.float64_t, ndim=1] wval,
                         method):
        if method == "mip" or method == -1:
            self.merged = -1
        elif method == "integrate" or method == 1:
            self.merged = 1
        cdef int curpos = 0
        self.num_cells = wval.shape[0]
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                curpos = self.unfill_buffer(self.root_nodes[i][j], curpos,
                                 refined, values, wval)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def tobuffer(self):
        cdef int total = self.num_cells
        # We now have four buffers:
        # Refined or not (total,) int32
        # Values in each node (total, nvals) float64
        # Weight values in each node (total,) float64
        cdef np.ndarray[np.int32_t, ndim=1] refined
        refined = np.zeros(total, dtype='int32')
        cdef np.ndarray[np.float64_t, ndim=2] values
        values = np.zeros((total, self.nvals), dtype='float64')
        cdef np.ndarray[np.float64_t, ndim=1] wval
        wval = np.zeros(total, dtype='float64')
        cdef int curpos = 0
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                curpos = self.fill_buffer(self.root_nodes[i][j], curpos,
                                 refined, values, wval)
        return (refined, values, wval)

    def get_args(self):
        return (self.top_grid_dims[0], self.top_grid_dims[1], self.nvals)

    cdef int add_to_position(self,
                 int level, np.int64_t pos[2],
                 np.float64_t *val,
                 np.float64_t weight_val, int skip = 0):
        cdef int i, j, L
        cdef QuadTreeNode *node
        node = self.find_on_root_level(pos, level)
        if node == NULL:
            return -1
        if level > self.max_level:
            self.max_level = level
        for L in range(level):
            if node.children[0][0] == NULL:
                QTN_refine(node, self.nvals)
                self.num_cells += 4
            # Maybe we should use bitwise operators?
            i = (pos[0] >> (level - L - 1)) & 1
            j = (pos[1] >> (level - L - 1)) & 1
            node = node.children[i][j]
        if skip == 1: return 0
        self.combine(node, val, weight_val, self.nvals)
        return 0

    @cython.cdivision(True)
    cdef QuadTreeNode *find_on_root_level(self, np.int64_t pos[2], int level):
        # We need this because the root level won't just have four children
        # So we find on the root level, then we traverse the tree.
        cdef np.int64_t i, j
        i = pos[0] >> level
        j = pos[1] >> level
        if i >= self.top_grid_dims[0] or i < 0 or \
           j >= self.top_grid_dims[1] or j < 0:
            self.last_dims[0] = i
            self.last_dims[1] = j
            return NULL
        return self.root_nodes[i][j]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def add_array_to_tree(self, int level, np.ndarray[np.int64_t, ndim=1] pxs,
            np.ndarray[np.int64_t, ndim=1] pys,
            np.ndarray[np.float64_t, ndim=2] pvals,
            np.ndarray[np.float64_t, ndim=1] pweight_vals,
            int skip = 0):
        cdef int p
        cdef np.float64_t *vals
        cdef np.float64_t *data = <np.float64_t *> pvals.data
        cdef np.int64_t pos[2]
        for p in range(pxs.shape[0]):
            vals = data + self.nvals*p
            pos[0] = pxs[p]
            pos[1] = pys[p]
            self.add_to_position(level, pos, vals, pweight_vals[p], skip)
        return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def add_chunk_to_tree(self,
            np.ndarray[np.int64_t, ndim=1] pxs,
            np.ndarray[np.int64_t, ndim=1] pys,
            np.ndarray[np.int64_t, ndim=1] level,
            np.ndarray[np.float64_t, ndim=2] pvals,
            np.ndarray[np.float64_t, ndim=1] pweight_vals):
        cdef int ps = pxs.shape[0]
        cdef int p, rv
        cdef np.float64_t *vals
        cdef np.float64_t *data = <np.float64_t *> pvals.data
        cdef np.int64_t pos[2]
        for p in range(ps):
            vals = data + self.nvals*p
            pos[0] = pxs[p]
            pos[1] = pys[p]
            rv = self.add_to_position(level[p], pos, vals, pweight_vals[p])
            if rv == -1:
                raise YTIntDomainOverflow(
                    (self.last_dims[0], self.last_dims[1]),
                    (self.top_grid_dims[0], self.top_grid_dims[1]))
        return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def initialize_chunk(self,
            np.ndarray[np.int64_t, ndim=1] pxs,
            np.ndarray[np.int64_t, ndim=1] pys,
            np.ndarray[np.int64_t, ndim=1] level):
        cdef int num = pxs.shape[0]
        cdef int p, rv
        cdef np.int64_t pos[2]
        for p in range(num):
            pos[0] = pxs[p]
            pos[1] = pys[p]
            rv = self.add_to_position(level[p], pos, NULL, 0.0, 1)
            if rv == -1:
                raise YTIntDomainOverflow(
                    (self.last_dims[0], self.last_dims[1]),
                    (self.top_grid_dims[0], self.top_grid_dims[1]))
        return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def get_all(self, int count_only = 0, int method = 1):
        cdef int i, j, vi
        cdef int total = 0
        self.merged = method
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                total += self.count(self.root_nodes[i][j])
        if count_only: return total
        # Allocate our array
        cdef np.ndarray[np.float64_t, ndim=1] opx
        cdef np.ndarray[np.float64_t, ndim=1] opy
        cdef np.ndarray[np.float64_t, ndim=1] opdx
        cdef np.ndarray[np.float64_t, ndim=1] opdy
        cdef np.ndarray[np.float64_t, ndim=2] nvals
        cdef np.ndarray[np.float64_t, ndim=1] nwvals
        opx = np.zeros(total, dtype='float64')
        opy = np.zeros(total, dtype='float64')
        opdx = np.zeros(total, dtype='float64')
        opdy = np.zeros(total, dtype='float64')
        nvals = np.zeros( (total, self.nvals), dtype='float64')
        nwvals = np.zeros( total, dtype='float64')
        cdef np.int64_t curpos = 0
        cdef np.float64_t *px = <np.float64_t *> opx.data
        cdef np.float64_t *py = <np.float64_t *> opy.data
        cdef np.float64_t *pdx = <np.float64_t *> opdx.data
        cdef np.float64_t *pdy = <np.float64_t *> opdy.data
        cdef np.float64_t *vdata = <np.float64_t *> nvals.data
        cdef np.float64_t *wdata = <np.float64_t *> nwvals.data
        cdef np.float64_t wtoadd
        cdef np.float64_t *vtoadd = <np.float64_t *> malloc(
                sizeof(np.float64_t)*self.nvals)
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                for vi in range(self.nvals): vtoadd[vi] = 0.0
                wtoadd = 0.0
                curpos += self.fill(self.root_nodes[i][j],
                    curpos, px, py, pdx, pdy, vdata, wdata, vtoadd, wtoadd, 0)
        free(vtoadd)
        return opx, opy, opdx, opdy, nvals, nwvals

    cdef int count(self, QuadTreeNode *node):
        cdef int i, j
        cdef int count = 0
        if node.children[0][0] == NULL: return 1
        for i in range(2):
            for j in range(2):
                count += self.count(node.children[i][j])
        return count

    @cython.cdivision(True)
    cdef int fill(self, QuadTreeNode *node,
                        np.int64_t curpos,
                        np.float64_t *px,
                        np.float64_t *py,
                        np.float64_t *pdx,
                        np.float64_t *pdy,
                        np.float64_t *vdata,
                        np.float64_t *wdata,
                        np.float64_t *vtoadd,
                        np.float64_t wtoadd,
                        np.int64_t level):
        cdef int i, j, n
        cdef np.float64_t *vorig
        vorig = <np.float64_t *> malloc(sizeof(np.float64_t) * self.nvals)
        if node.children[0][0] == NULL:
            if self.merged == -1:
                for i in range(self.nvals):
                    vdata[self.nvals * curpos + i] = fmax(node.val[i], vtoadd[i])
                wdata[curpos] = 1.0
            else:
                for i in range(self.nvals):
                    vdata[self.nvals * curpos + i] = node.val[i] + vtoadd[i]
                wdata[curpos] = node.weight_val + wtoadd
            pdx[curpos] = 1.0 / (self.top_grid_dims[0]*2**level)
            pdy[curpos] = 1.0 / (self.top_grid_dims[1]*2**level)
            px[curpos] = (0.5 + node.pos[0]) * pdx[curpos]
            py[curpos] = (0.5 + node.pos[1]) * pdy[curpos]
            pdx[curpos] /= 2.0
            pdy[curpos] /= 2.0
            return 1
        cdef np.int64_t added = 0
        if self.merged == 1:
            for i in range(self.nvals):
                vorig[i] = vtoadd[i]
                vtoadd[i] += node.val[i]
            wtoadd += node.weight_val
        elif self.merged == -1:
            for i in range(self.nvals):
                vtoadd[i] = node.val[i]
        for i in range(2):
            for j in range(2):
                if self.merged == -1:
                    for n in range(self.nvals):
                        vtoadd[n] = node.val[n]
                added += self.fill(node.children[i][j],
                        curpos + added, px, py, pdx, pdy, vdata, wdata,
                        vtoadd, wtoadd, level + 1)
        if self.merged == 1:
            for i in range(self.nvals):
                vtoadd[i] = vorig[i]
            wtoadd -= node.weight_val
        free(vorig)
        return added

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_image(self, np.ndarray[np.float64_t, ndim=2] buffer, _bounds,
                   int val_index = 0, int weighted = 0):
        cdef np.float64_t pos[2]
        cdef np.float64_t dds[2]
        cdef int nn[2]
        cdef int i, j
        cdef np.float64_t bounds[4]
        cdef np.float64_t opos[4]
        cdef np.float64_t weight = 0.0, value = 0.0
        cdef np.float64_t *wval = NULL
        if weighted == 1:
            wval = &weight
        for i in range(4):
            bounds[i] = _bounds[i]
        for i in range(2):
            nn[i] = buffer.shape[i]
            dds[i] = (bounds[i*2 + 1] - bounds[i*2])/nn[i]
        pos[0] = bounds[0]
        opos[0] = opos[1] = pos[0] + dds[0]
        for i in range(nn[0]):
            pos[1] = bounds[2]
            opos[2] = opos[3] = pos[1] + dds[1]
            for j in range(nn[1]):
                # We start at level zero.  In the future we could optimize by
                # retaining oct information from previous cells.
                if not ((opos[0] <= pos[0] <= opos[1]) and
                        (opos[2] <= pos[1] <= opos[3])):
                    value = self.find_value_at_pos(pos, val_index,
                                        opos, wval)
                buffer[i,j] = value
                if weighted == 1:
                    buffer[i,j] /= weight
                pos[1] += dds[1]
            pos[0] += dds[0]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.float64_t find_value_at_pos(self, np.float64_t pos[2],
                                         int val_index, np.float64_t opos[4],
                                         np.float64_t *wval = NULL):
        cdef int i
        cdef np.int64_t ind[2]
        cdef np.float64_t cp[2]
        cdef np.float64_t dds[2]
        cdef np.float64_t value = 0.0
        cdef np.float64_t weight = 0.0
        cdef QuadTreeNode *cur
        for i in range(2):
            ind[i] = <np.int64_t> (pos[i]/self.dds[i])
            cp[i] = (ind[i] + 0.5) * self.dds[i]
            dds[i] = self.dds[i]
        cur = self.root_nodes[ind[0]][ind[1]]
        value += cur.val[val_index]
        weight += cur.weight_val
        while cur.children[0][0] != NULL:
            for i in range(2):
                # Note that below offset by half a dx for center, after
                # updating to the next level
                dds[i] = dds[i] * 0.5
                if cp[i] >= pos[i]:
                    ind[i] = 0
                    cp[i] -= dds[i] * 0.5
                else:
                    ind[i] = 1
                    cp[i] += dds[i] * 0.5
            cur = cur.children[ind[0]][ind[1]]
            value += cur.val[val_index]
            weight += cur.weight_val
        opos[0] = cp[0] - dds[0] * 0.5
        opos[1] = cp[0] + dds[0] * 0.5
        opos[2] = cp[1] - dds[1] * 0.5
        opos[3] = cp[1] + dds[1] * 0.5
        if wval != NULL:
            wval[0] = weight
        return value

    def __dealloc__(self):
        cdef int i, j
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                QTN_free(self.root_nodes[i][j])
            free(self.root_nodes[i])
        free(self.root_nodes)

cdef void QTN_merge_nodes(QuadTreeNode *n1, QuadTreeNode *n2, int nvals,
                          QTN_combine *func):
    # We have four choices when merging nodes.
    # 1. If both nodes have no refinement, then we add values of n2 to n1.
    # 2. If both have refinement, we call QTN_merge_nodes on all four children.
    # 3. If n2 has refinement and n1 does not, we detach n2's children and
    #    attach them to n1.
    # 4. If n1 has refinement and n2 does not, we add the value of n2 to n1.
    cdef int i, j

    func(n1, n2.val, n2.weight_val, nvals)
    if n1.children[0][0] == n2.children[0][0] == NULL:
        pass
    elif n1.children[0][0] != NULL and n2.children[0][0] != NULL:
        for i in range(2):
            for j in range(2):
                QTN_merge_nodes(n1.children[i][j], n2.children[i][j], nvals, func)
    elif n1.children[0][0] == NULL and n2.children[0][0] != NULL:
        for i in range(2):
            for j in range(2):
                n1.children[i][j] = n2.children[i][j]
                n2.children[i][j] = NULL
    elif n1.children[0][0] != NULL and n2.children[0][0] == NULL:
        pass
    else:
        raise RuntimeError

def merge_quadtrees(QuadTree qt1, QuadTree qt2, method = 1):
    cdef int i, j
    qt1.num_cells = 0
    cdef QTN_combine *func
    if method == 1:
        qt1.merged = 1
        func = QTN_add_value
    elif method == -1:
        qt1.merged = -1
        func = QTN_max_value
    else:
        raise NotImplementedError
    if qt1.merged != 0 or qt2.merged != 0:
        assert(qt1.merged == qt2.merged)
    for i in range(qt1.top_grid_dims[0]):
        for j in range(qt1.top_grid_dims[1]):
            QTN_merge_nodes(qt1.root_nodes[i][j],
                            qt2.root_nodes[i][j],
                            qt1.nvals, func)
            qt1.num_cells += qt1.count_total_cells(
                                qt1.root_nodes[i][j])
