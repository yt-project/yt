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

cdef inline np.float64_t f64max(np.float64_t f0, np.float64_t f1):
    if f0 > f1: return f0
    return f1

cdef struct OctreeNode:
    np.float64_t *val
    np.float64_t weight_val
    np.int64_t pos[3]
    int level
    int nvals
    OctreeNode *children[2][2][2]
    OctreeNode *parent

cdef void OTN_add_value(OctreeNode *self,
        np.float64_t *val, np.float64_t weight_val):
    cdef int i
    for i in range(self.nvals):
        self.val[i] += val[i]
    self.weight_val += weight_val

cdef void OTN_refine(OctreeNode *self, int incremental = 0):
    cdef int i, j, k, i1, j1
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
                            self.level + 1, self)
    if incremental: return
    for i in range(self.nvals): self.val[i] = 0.0
    self.weight_val = 0.0

cdef int OTN_same(OctreeNode *node1, OctreeNode *node2):
    # Returns 1 if node1 == node2; 0 otherwise.
    # If all of these conditions are met, they are the same node.
    if node1.pos[0] == node2.pos[0] and \
       node1.pos[1] == node2.pos[1] and \
       node1.pos[2] == node2.pos[2] and \
       node1.level == node2.level: return 1
    return 0

cdef int OTN_contained(OctreeNode *node1, OctreeNode *node2):
    # Returns 1 if node2 contains node1; and 0 otherwise.
    # node1.level > node2.level.
    cdef OctreeNode *parent_node
    parent_node = node1.parent
    while parent_node is not NULL:
        if parent_node is node2: return 1
        parent_node = parent_node.parent
    # If we've gotten this far, the two nodes are not related.
    return 0

cdef OctreeNode *OTN_initialize(np.int64_t pos[3], int nvals,
                        np.float64_t *val, np.float64_t weight_val,
                        int level, OctreeNode *parent):
    cdef OctreeNode *node
    cdef int i, j, k
    node = <OctreeNode *> malloc(sizeof(OctreeNode))
    node.pos[0] = pos[0]
    node.pos[1] = pos[1]
    node.pos[2] = pos[2]
    node.nvals = nvals
    node.parent = parent
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
    cdef int i, j, k
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
    cdef int incremental
    cdef np.float64_t opening_angle
    cdef int count1, count2

    def __cinit__(self, np.ndarray[np.int64_t, ndim=1] top_grid_dims,
                  int nvals, int incremental = False):
        cdef int i, j, k
        self.incremental = incremental
        cdef OctreeNode *node
        cdef np.int64_t pos[3]
        cdef np.float64_t *vals = <np.float64_t *> alloca(
                sizeof(np.float64_t)*nvals)
        cdef np.float64_t weight_val = 0.0
        self.nvals = nvals
        count1 = 0
        count2 = 0
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
                        pos, nvals, vals, weight_val, 0, NULL)

    cdef void add_to_position(self,
                 int level, np.int64_t pos[3],
                 np.float64_t *val,
                 np.float64_t weight_val):
        cdef int i, j, k, L
        cdef OctreeNode *node
        node = self.find_on_root_level(pos, level)
        cdef np.int64_t fac
        for L in range(level):
            if self.incremental:
                OTN_add_value(node, val, weight_val)
            if node.children[0][0][0] == NULL:
                OTN_refine(node, self.incremental)
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
        cdef np.int64_t i, j, k
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
        cdef int i, j, k
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
        npos = np.zeros( (total, 3), dtype='int64')
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
        cdef int i, j, k
        # We only really return a non-zero, calculated value if we are at the
        # level in question.
        if node.level == level:
            if self.incremental: return 1
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
        cdef int i, j, k
        if node.level == level:
            if node.children[0][0][0] != NULL and not self.incremental:
                return 0
            for i in range(self.nvals):
                vdata[self.nvals * curpos + i] = node.val[i]
            wdata[curpos] = node.weight_val
            pdata[curpos * 3] = node.pos[0]
            pdata[curpos * 3 + 1] = node.pos[1]
            pdata[curpos * 3 + 2] = node.pos[2]
            return 1
        if node.children[0][0][0] == NULL: return 0
        cdef np.int64_t added = 0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    added += self.fill_from_level(node.children[i][j][k],
                            level, curpos + added, pdata, vdata, wdata)
        return added

    cdef np.float64_t fbe_node_separation(self, OctreeNode *node1, OctreeNode *node2):
        # Find the distance between the two nodes. To match FindBindingEnergy
        # in data_point_utilities.c, we'll do this in code [0,1) units.
        # We'll also assume 
        cdef np.float64_t fn1, fn2, p1, p2, dist
        cdef np.int64_t n1, n2
        cdef int i
        dist = 0.0
        for i in range(3):
            n1 = self.po2[node1.level] * self.top_grid_dims[i]
            n2 = self.po2[node2.level] * self.top_grid_dims[i]
            fn1 = <np.float64_t> n1
            fn2 = <np.float64_t> n2
            # The added term is to re-cell center the data.
            p1 = (<np.float64_t> node1.pos[i]) / fn1 + 1. / fn1 / 2.
            p2 = (<np.float64_t> node2.pos[i]) / fn2 + 1. / fn2 / 2.
            dist += (p1 - p2) * (p1 - p2)
        dist = np.sqrt(dist)
        return dist
    
    cdef np.float64_t fbe_opening_angle(self, OctreeNode *node1,
            OctreeNode *node2):
        # Calculate the opening angle of node2 upon the center of node1.
        # In order to keep things simple, we will not assume symmetry in all
        # three directions of the octree, and we'll use the largest dimension
        # if the tree is not symmetric. This is not strictly the opening angle
        # the purest sense, but it's slightly more accurate, so it's OK.
        # This is done in code units to match the distance calculation.
        cdef np.float64_t d2, dx2, dist
        cdef np.int64_t n2
        cdef int i
        d2 = 0.0
        if OTN_same(node1, node2): return 100000.0 # Just some large number.
        if self.top_grid_dims[1] == self.top_grid_dims[0] and \
                self.top_grid_dims[2] == self.top_grid_dims[0]:
            # Symmetric
            n2 = self.po2[node2.level] * self.top_grid_dims[0]
            d2 = 1. / (<np.float64_t> n2)
        else:
            # Not symmetric
            for i in range(3):
                n2 = self.po2[node2.level] * self.top_grid_dims[i]
                dx2 = 1. / (<np.float64_t> n2)
                d2 = f64max(d2, dx2)
        # Now calculate the opening angle.
        dist = self.fbe_node_separation(node1, node2)
        return d2 / dist

    cdef np.float64_t fbe_potential_of_remote_nodes(self, OctreeNode *node1,
            int top_i, int top_j, int top_k):
        # Given a childless node "node1", calculate the potential for it from
        # all the other nodes using the treecode method.
        cdef int i, j, k
        cdef np.float64_t potential
        potential = 0.0
        # We *do* want to walk over the root_node that node1 is in to look for
        # level>0 nodes close to node1, but none of
        # the previously-walked root_nodes.
        for i in range(top_i, self.top_grid_dims[0]):
            for j in range(top_j, self.top_grid_dims[1]):
                for k in range(top_k, self.top_grid_dims[2]):
                    if self.root_nodes[i][j][k].val[0] == 0: continue
                    potential += self.fbe_iterate_remote_nodes(node1,
                        self.root_nodes[i][j][k])
        return potential

    cdef np.float64_t fbe_iterate_remote_nodes(self, OctreeNode *node1,
            OctreeNode *node2):
        # node1 never changes.
        # node2 is the iterated-upon remote node.
        cdef int i, j, k, contained
        cdef np.float64_t potential, dist, angle
        potential = 0.0
        # Do nothing with node2 when it is the same as node1. No potential
        # calculation, and no digging deeper.
        if OTN_same(node1, node2): return 0.0
        # Is node1 contained inside node2?
        contained = OTN_contained(node1, node2)
        # If we have a childless node2 we can skip everything below and
        # calculate the potential, as long as node2 does not contain node1.
        if node2.children[0][0][0] == NULL and not contained:
            dist = self.fbe_node_separation(node1, node2)
            self.count2 += 1
            #print '%1.10e %1.10e %1.10e %1.10e xx' % (node1.val[0], node2.val[0],
            #dist, node1.val[0] * node2.val[0] / dist)
            return node1.val[0] * node2.val[0] / dist
        # Now we apply the opening angle test. If the opening angle is small
        # enough, we use this node for the potential and dig no deeper.
        angle = self.fbe_opening_angle(node1, node2)
        if angle < self.opening_angle and not contained:
            dist = self.fbe_node_separation(node1, node2)
            return node1.val[0] * node2.val[0] / dist
        # If the above is not satisfied, we must dig deeper!
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if node2.children[i][j][k].val[0] == 0.0: continue
                    potential += self.fbe_iterate_remote_nodes(node1,
                        node2.children[i][j][k])
        return potential

    cdef np.float64_t fbe_iterate_children(self, OctreeNode *node, int top_i,
            int top_j, int top_k):
        # Recursively iterate over child nodes until we get a childless node.
        cdef int i, j, k
        cdef np.float64_t potential
        potential = 0.0
        # We have a childless node. Time to iterate over every other
        # node using the treecode method.
        if node.children[0][0][0] is NULL:
            self.count1 += 1
            potential = self.fbe_potential_of_remote_nodes(node,
                top_i, top_j, top_k)
            return potential
        # If the node has children, we need to walk all of them returning
        # the potential for each.
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if node.children[i][j][k].val[0] == 0.0: continue
                    potential += self.fbe_iterate_children(node.children[i][j][k],
                        top_i, top_j, top_k)
        return potential

    def find_binding_energy(self, truncate, kinetic, opening_angle = 1.0):
        r"""Find the binding energy of an ensemble of data points using the
        treecode method.
        
        Note: The first entry of the vals array MUST be Mass.
        """
        # Here are the order of events:
        # 1. We loop over all of the root_nodes, below.
        # 2. In fbe_iterate_children, each of these nodes is iterated until
        #    we reach a node without any children.
        # 3. Next, starting in fbe_potential_of_remote_nodes we again loop over
        #    root_nodes, not double-counting, calling fbe_iterate_remote_nodes.
        # 4. In fbe_iterate_remote_nodes, if we have a childless node, we
        #    calculate the potential between the two nodes and return. Or, if
        #    we have a node with a small enough opening angle, we return the
        #    potential. If neither of these are done, we call
        #    fbe_iterate_remote_nodes on the children nodes.
        # 5. All of this returns a total, non-double-counted potential.
        cdef int i, j, k
        cdef np.float64_t potential
        potential = 0.0
        self.opening_angle = opening_angle
        # The first part of the loop goes over all of the root level nodes.
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                for k in range(self.top_grid_dims[2]):
                    if self.root_nodes[i][j][k].val[0] == 0.0: continue
                    #print i,j,k
                    potential += self.fbe_iterate_children(self.root_nodes[i][j][k],
                        i, j, k)
                    #potential = kinetic * 2
                    if truncate and potential > kinetic: break
                if truncate and potential > kinetic: break
            if truncate and potential > kinetic:
                print "Truncating!"
                break
        print 'count1', self.count1
        print 'count2', self.count2
        return potential

    def __dealloc__(self):
        cdef int i, j, k
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                for k in range(self.top_grid_dims[2]):
                    OTN_free(self.root_nodes[i][j][k])
                free(self.root_nodes[i][j])
            free(self.root_nodes[i])
        free(self.root_nodes)

