"""
A refine-by-two AMR-specific octree



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
# Double up here for def'd functions
cimport numpy as cnp
cimport cython

from yt.utilities.lib.fp_utils cimport imax, fmax, imin, fmin, iclip, fclip
from libc.stdlib cimport malloc, free, abs

import sys, time

cdef extern from "platform_dep.h":
    # NOTE that size_t might not be int
    void *alloca(int)

cdef extern from "math.h":
    double sqrt(double x)

cdef inline np.float64_t f64max(np.float64_t f0, np.float64_t f1):
    if f0 > f1: return f0
    return f1

cdef struct OctreeNode:
    np.float64_t *val
    np.float64_t weight_val
    np.int64_t pos[3]
    int level
    int nvals
    int max_level # The maximum level under this node with mass.
    OctreeNode *children[2][2][2]
    OctreeNode *parent
    OctreeNode *next
    OctreeNode *up_next

cdef void OTN_add_value(OctreeNode *self,
        np.float64_t *val, np.float64_t weight_val, int level, int treecode):
    cdef int i
    for i in range(self.nvals):
        self.val[i] += val[i]
    self.weight_val += weight_val
    if treecode and val[0] > 0.:
        self.max_level = imax(self.max_level, level)

cdef void OTN_refine(OctreeNode *self, int incremental = 0):
    cdef int i, j, k
    cdef np.int64_t npos[3]
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
                            self.level + 1, self, incremental)
    if incremental: return
    for i in range(self.nvals): self.val[i] = 0.0
    self.weight_val = 0.0

cdef OctreeNode *OTN_initialize(np.int64_t pos[3], int nvals,
                        np.float64_t *val, np.float64_t weight_val,
                        int level, OctreeNode *parent, int incremental = 0):
    cdef OctreeNode *node
    cdef int i, j, k
    node = <OctreeNode *> malloc(sizeof(OctreeNode))
    node.pos[0] = pos[0]
    node.pos[1] = pos[1]
    node.pos[2] = pos[2]
    node.nvals = nvals
    node.parent = parent
    node.next = NULL
    node.up_next = NULL
    node.max_level = 0
    node.val = <np.float64_t *> malloc(
                nvals * sizeof(np.float64_t))
    if incremental:
        for i in range(nvals):
            node.val[i] = 0.
        node.weight_val = 0.
    else:
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
    # Below is for the treecode.
    cdef np.float64_t opening_angle
    # We'll store dist here so it doesn't have to be calculated twice.
    cdef np.float64_t dist
    cdef np.float64_t root_dx[3]
    cdef OctreeNode *last_node

    def __cinit__(self, np.ndarray[np.int64_t, ndim=1] top_grid_dims,
                  int nvals, int incremental = False):
        cdef int i, j, k
        self.incremental = incremental
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
                        pos, nvals, vals, weight_val, 0, NULL)

    cdef void add_to_position(self,
                 int level, np.int64_t pos[3],
                 np.float64_t *val,
                 np.float64_t weight_val, treecode):
        cdef int i, j, k, L
        cdef OctreeNode *node
        node = self.find_on_root_level(pos, level)
        cdef np.int64_t fac
        for L in range(level):
            if self.incremental:
                OTN_add_value(node, val, weight_val, level, treecode)
            if node.children[0][0][0] == NULL:
                OTN_refine(node, self.incremental)
            # Maybe we should use bitwise operators?
            fac = self.po2[level - L - 1]
            i = (pos[0] >= fac*(2*node.pos[0]+1))
            j = (pos[1] >= fac*(2*node.pos[1]+1))
            k = (pos[2] >= fac*(2*node.pos[2]+1))
            node = node.children[i][j][k]
        OTN_add_value(node, val, weight_val, level, treecode)
            
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
            np.ndarray[np.float64_t, ndim=1] pweight_vals,
            int treecode = 0):
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
            self.add_to_position(level, pos, vals, pweight_vals[p], treecode)

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

    @cython.cdivision(True)
    cdef np.float64_t fbe_node_separation(self, OctreeNode *node1, OctreeNode *node2):
        # Find the distance between the two nodes.
        cdef np.float64_t dx1, dx2, p1, p2, dist
        cdef int i
        dist = 0.0
        for i in range(3):
            # Discover the appropriate dx for each node/dim.
            dx1 = self.root_dx[i] / (<np.float64_t> self.po2[node1.level])
            dx2 = self.root_dx[i] / (<np.float64_t> self.po2[node2.level])
            # The added term is to re-cell center the data.
            p1 = (<np.float64_t> node1.pos[i]) * dx1 + dx1/2.
            p2 = (<np.float64_t> node2.pos[i]) * dx2 + dx2/2.
            dist += (p1 - p2) * (p1 - p2)
        dist = sqrt(dist)
        return dist
    
    @cython.cdivision(True)
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
        if node1 is node2: return 100000.0 # Just some large number.
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
        self.dist = dist
        return d2 / dist

    cdef void set_next(self, OctreeNode *node, int treecode):
        # This sets up the linked list, pointing node.next to the next node
        # in the iteration order.
        cdef int i, j, k
        if treecode and node.val[0] is not 0.:
            # In a treecode, we only make a new next link if this node has mass.
            self.last_node.next = node
            self.last_node = node
        elif treecode and node.val[0] is 0.:
            # No mass means it's children have no mass, no need to dig an
            # further.
            return
        else:
            # We're not doing the treecode, but we still want a linked list,
            # we don't care about val[0] necessarily.
            self.last_node.next = node
            self.last_node = node
        if node.children[0][0][0] is NULL: return
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    self.set_next(node.children[i][j][k], treecode)
        return

    cdef void set_up_next(self, OctreeNode *node):
        # This sets up a second linked list, pointing node.up_next to the next
        # node in the list that is at the same or lower (coarser) level than
        # this node. This is useful in the treecode for skipping over nodes
        # that don't need to be inspected.
        cdef OctreeNode *initial_next
        cdef OctreeNode *temp_next
        initial_next = node.next
        temp_next = node.next
        if node.next is NULL: return
        while temp_next.level > node.level:
            temp_next = temp_next.next
            if temp_next is NULL: break
        node.up_next = temp_next
        self.set_up_next(initial_next)

    def finalize(self, int treecode = 0):
        # Set up the linked list for the nodes.
        # Set treecode = 1 if nodes with no mass are to be skipped in the
        # list.
        cdef int i, j, k, sum, top_grid_total, ii, jj, kk
        self.last_node = self.root_nodes[0][0][0]
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                for k in range(self.top_grid_dims[2]):
                    self.set_next(self.root_nodes[i][j][k], treecode)
        # Now we want to link to the next node in the list that is
        # on a level the same or lower (coarser) than us. This will be used
        # during a treecode search so we can skip higher-level (finer) nodes.
        sum = 1
        top_grid_total = self.top_grid_dims[0] * self.top_grid_dims[1] * \
            self.top_grid_dims[2]
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                for k in range(self.top_grid_dims[2]):
                    self.set_up_next(self.root_nodes[i][j][k])
                    # Point the root_nodes.up_next to the next root_node in the
                    # list, except for the last one which stays pointing to NULL.
                    if sum < top_grid_total - 1:
                        ii = i
                        jj = j
                        kk = (k + 1) % self.top_grid_dims[2]
                        if kk < k:
                            jj = (j + 1) % self.top_grid_dims[1]
                            if jj < j:
                                ii = (i + 1) % self.top_grid_dims[0]
                        self.root_nodes[i][j][k].up_next = \
                            self.root_nodes[ii][jj][kk]
                    sum += 1

    @cython.cdivision(True)
    cdef np.float64_t fbe_main(self, np.float64_t potential, int truncate,
            np.float64_t kinetic):
        # The work is done here. Starting at the top of the linked list of
        # nodes, 
        cdef np.float64_t angle, dist
        cdef OctreeNode *this_node
        cdef OctreeNode *pair_node
        this_node = self.root_nodes[0][0][0]
        while this_node is not NULL:
            # Iterate down the list to a node that either has no children and
            # is at the max_level of the tree, or to a node where
            # all of its children are massless. The second case is when data
            # from a level that isn't the deepest has been added to the tree.
            while this_node.max_level is not this_node.level:
                this_node = this_node.next
                # In case we reach the end of the list...
                if this_node is NULL: break
            if this_node is NULL: break
            if truncate and potential > kinetic:
                print 'Truncating...'
                break
            pair_node = this_node.next
            while pair_node is not NULL:
                # If pair_node is massless, we can jump to up_next, because
                # nothing pair_node contains will have mass either.
                # I think that this should primarily happen for root_nodes
                # created for padding to make the region cubical.
                if pair_node.val[0] is 0.0:
                    pair_node = pair_node.up_next
                    continue
                # If pair_node is a childless node, or is a coarser node with
                # no children, we can calculate the pot
                # right now, and get a new pair_node.
                if pair_node.max_level is pair_node.level:
                    dist = self.fbe_node_separation(this_node, pair_node)
                    potential += this_node.val[0] * pair_node.val[0] / dist
                    if truncate and potential > kinetic: break
                    pair_node = pair_node.next
                    continue
                # Next, if the opening angle to pair_node is small enough,
                # calculate the potential and get a new pair_node using
                # up_next because we don't need to look at pair_node's children.
                angle = self.fbe_opening_angle(this_node, pair_node)
                if angle < self.opening_angle:
                    # self.dist was just set in fbe_opening_angle, so we
                    # can use it here without re-calculating it for these two
                    # nodes.
                    potential += this_node.val[0] * pair_node.val[0] / self.dist
                    if truncate and potential > kinetic: break
                    # We can skip all the nodes that are contained within 
                    # pair_node, saving time walking the linked list.
                    pair_node = pair_node.up_next
                # If we've gotten this far, pair_node has children, but it's
                # too coarse, so we simply dig deeper using .next.
                else:
                    pair_node = pair_node.next
            # We've exhausted the pair_nodes.
            # Now we find a new this_node in the list, and do more searches
            # over pair_node.
            this_node = this_node.next
        return potential

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def find_binding_energy(self, int truncate, np.float64_t kinetic,
        np.ndarray[np.float64_t, ndim=1] root_dx, float opening_angle = 1.0):
        r"""Find the binding energy of an ensemble of data points using the
        treecode method.
        
        Note: The first entry of the vals array MUST be Mass. Any other
        values will be ignored, including the weight array.
        """
        # The real work is done in fbe_main(), this just sets things up
        # and returns the potential.
        cdef int i
        cdef np.float64_t potential
        potential = 0.0
        self.opening_angle = opening_angle
        for i in range(3):
            self.root_dx[i] = root_dx[i]
        potential = self.fbe_main(potential, truncate, kinetic)
        return potential

    cdef int node_ID(self, OctreeNode *node):
        # Returns an unique ID for this node based on its position and level.
        cdef int ID, i, offset, root
        cdef np.int64_t this_grid_dims[3]
        offset = 0
        root = 1
        for i in range(3):
            root *= self.top_grid_dims[i]
            this_grid_dims[i] = self.top_grid_dims[i] * 2**node.level
        for i in range(node.level):
            offset += root * 2**(3 * i)
        ID = offset + (node.pos[0] + this_grid_dims[0] * (node.pos[1] + \
            this_grid_dims[1] * node.pos[2]))
        return ID

    cdef int node_ID_on_level(self, OctreeNode *node):
        # Returns the node ID on node.level for this node.
        cdef int ID, i
        cdef np.int64_t this_grid_dims[3]
        for i in range(3):
            this_grid_dims[i] = self.top_grid_dims[i] * 2**node.level
        ID = node.pos[0] + this_grid_dims[0] * (node.pos[1] + \
            this_grid_dims[1] * node.pos[2])
        return ID

    cdef void print_node_info(self, OctreeNode *node):
        cdef int i, j, k
        line = "%d\t" % self.node_ID(node)
        if node.next is not NULL:
            line += "%d\t" % self.node_ID(node.next)
        else: line += "-1\t"
        if node.up_next is not NULL:
            line += "%d\t" % self.node_ID(node.up_next)
        else: line += "-1\t"
        line += "%d\t%d\t%d\t%d\t" % (node.level,node.pos[0],node.pos[1],node.pos[2])
        for i in range(node.nvals):
            line += "%1.5e\t" % node.val[i]
        line += "%f\t" % node.weight_val
        line += "%s\t%s\t" % (node.children[0][0][0] is not NULL, node.parent is not NULL)
        if node.children[0][0][0] is not NULL:
            nline = ""
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        nline += "%d," % self.node_ID(node.children[i][j][k])
            line += nline
        print line
        return

    cdef void iterate_print_nodes(self, OctreeNode *node):
        cdef int i, j, k
        self.print_node_info(node)
        if node.children[0][0][0] is NULL:
            return
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    self.iterate_print_nodes(node.children[i][j][k])
        return

    def print_all_nodes(self):
        r"""
        Prints out information about all the nodes in the octree.
        
        Parameters
        ----------
        None.
        
        Examples
        --------
        >>> octree.print_all_nodes()
        (many lines of data)
        """
        cdef int i, j, k
        sys.stdout.flush()
        sys.stderr.flush()
        line = "ID\tnext\tup_n\tlevel\tx\ty\tz\t"
        for i in range(self.nvals):
            line += "val%d\t\t" % i
        line += "weight\t\tchild?\tparent?\tchildren"
        print line
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                for k in range(self.top_grid_dims[2]):
                    self.iterate_print_nodes(self.root_nodes[i][j][k])
        sys.stdout.flush()
        sys.stderr.flush()
        return

    def __dealloc__(self):
        cdef int i, j, k
        for i in range(self.top_grid_dims[0]):
            for j in range(self.top_grid_dims[1]):
                for k in range(self.top_grid_dims[2]):
                    OTN_free(self.root_nodes[i][j][k])
                free(self.root_nodes[i][j])
            free(self.root_nodes[i])
        free(self.root_nodes)

