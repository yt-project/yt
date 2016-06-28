"""
AMR kD-Tree Cython Tools



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
from libc.stdlib cimport malloc, free
from cython.view cimport array as cvarray

DEF Nch = 4

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef class Node:

    def __cinit__(self,
                  Node parent,
                  Node left,
                  Node right,
                  np.float64_t[:] left_edge,
                  np.float64_t[:] right_edge,
                  int grid,
                  np.int64_t node_id):
        self.dirty = False
        self.left = left
        self.right = right
        self.parent = parent
        cdef int i
        for i in range(3):
            self.left_edge[i] = left_edge[i]
            self.right_edge[i] = right_edge[i]
        self.grid = grid
        self.node_id = node_id
        self.split == NULL


    def print_me(self):
        print 'Node %i' % self.node_id
        print '\t le: %e %e %e' % (self.left_edge[0], self.left_edge[1],
                                   self.left_edge[2])
        print '\t re: %e %e %e' % (self.right_edge[0], self.right_edge[1],
                                   self.right_edge[2])
        print '\t grid: %i' % self.grid

    def get_split_dim(self):
        if self.split != NULL:
            return self.split.dim
        else:
            return -1

    def get_split_pos(self):
        if self.split != NULL:
            return self.split.pos
        else:
            return np.nan

    def set_left_edge(self, np.float64_t[:] left_edge):
        cdef int i
        for i in range(3):
            self.left_edge[i] = left_edge[i]

    def set_right_edge(self, np.float64_t[:] right_edge):
        cdef int i
        for i in range(3):
            self.right_edge[i] = right_edge[i]

    def create_split(self, dim, pos):
        split = <Split *> malloc(sizeof(Split))
        split.dim = dim
        split.pos = pos
        self.split = split

    def __dealloc__(self):
        if self.split != NULL: free(self.split)

    # Begin input of converted methods

    def get_left_edge(self):
        le = np.empty(3, dtype='float64')
        for i in range(3):
            le[i] = self.left_edge[i]
        return le

    def get_right_edge(self):
        re = np.empty(3, dtype='float64')
        for i in range(3):
            re[i] = self.right_edge[i]
        return re

    def set_dirty(self, bint state):
        cdef Node node
        for node in self.depth_traverse():
            node.dirty = state

    def kd_traverse(self, viewpoint=None):
        cdef Node node
        if viewpoint is None:
            for node in self.depth_traverse():
                if node._kd_is_leaf() == 1 and node.grid != -1:
                    yield node
        else:
            for node in self.viewpoint_traverse(viewpoint):
                if node._kd_is_leaf() == 1 and node.grid != -1:
                    yield node

    def add_pygrid(self,
                       np.float64_t[:] gle,
                       np.float64_t[:] gre,
                       int gid,
                       int rank,
                       int size):

        """
        The entire purpose of this function is to move everything from ndarrays
        to internal C pointers.
        """
        cdef np.float64_t[:] pgles = cvarray(format="d", shape=(3,), itemsize=sizeof(np.float64_t))
        cdef np.float64_t[:] pgres = cvarray(format="d", shape=(3,), itemsize=sizeof(np.float64_t))
        cdef int j
        for j in range(3):
            pgles[j] = gle[j]
            pgres[j] = gre[j]

        self.add_grid(pgles, pgres, gid, rank, size)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef add_grid(self,
                       np.float64_t[:] gle,
                       np.float64_t[:] gre,
                       int gid,
                       int rank,
                       int size):

        if not should_i_build(self, rank, size):
            return

        if self._kd_is_leaf() == 1:
            self.insert_grid(gle, gre, gid, rank, size)
        else:
            less_id = gle[self.split.dim] < self.split.pos
            if less_id:
                self.left.add_grid(gle, gre,
                         gid, rank, size)

            greater_id = gre[self.split.dim] > self.split.pos
            if greater_id:
                self.right.add_grid(gle, gre,
                         gid, rank, size)
        return



    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef insert_grid(self,
                    np.float64_t[:] gle,
                    np.float64_t[:] gre,
                    int grid_id,
                    int rank,
                    int size):
        if not should_i_build(self, rank, size):
            return

        # If we should continue to split based on parallelism, do so!
        if self.should_i_split(rank, size):
            self.geo_split(gle, gre, grid_id, rank, size)
            return

        cdef int contained = 1
        for i in range(3):
            if gle[i] > self.left_edge[i] or\
               gre[i] < self.right_edge[i]:
                contained *= 0

        if contained == 1:
            self.grid = grid_id
            assert(self.grid != -1)
            return

        # Split the grid
        cdef int check = self.split_grid(gle, gre, grid_id, rank, size)
        # If check is -1, then we have found a place where there are no choices.
        # Exit out and set the node to None.
        if check == -1:
            self.grid = -1
        return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef add_grids(self,
                        int ngrids,
                        np.float64_t[:,:] gles,
                        np.float64_t[:,:] gres,
                        np.int64_t[:] gids,
                        int rank,
                        int size):
        cdef int i, j, nless, ngreater, index
        cdef np.float64_t[:,:] less_gles, less_gres, greater_gles, greater_gres
        cdef np.int64_t[:] l_ids, g_ids
        if not should_i_build(self, rank, size):
            return

        if self._kd_is_leaf() == 1:
            self.insert_grids(ngrids, gles, gres, gids, rank, size)
            return

        less_ids = cvarray(format="q", shape=(ngrids,), itemsize=sizeof(np.int64_t))
        greater_ids = cvarray(format="q", shape=(ngrids,), itemsize=sizeof(np.int64_t))

        nless = 0
        ngreater = 0
        for i in range(ngrids):
            if gles[i, self.split.dim] < self.split.pos:
                less_ids[nless] = i
                nless += 1

            if gres[i, self.split.dim] > self.split.pos:
                greater_ids[ngreater] = i
                ngreater += 1

        #print 'nless: %i' % nless
        #print 'ngreater: %i' % ngreater

        if nless > 0:
            less_gles = cvarray(format="d", shape=(nless,3), itemsize=sizeof(np.float64_t))
            less_gres = cvarray(format="d", shape=(nless,3), itemsize=sizeof(np.float64_t))
            l_ids = cvarray(format="q", shape=(nless,), itemsize=sizeof(np.int64_t))

            for i in range(nless):
                index = less_ids[i]
                l_ids[i] = gids[index]
                for j in range(3):
                    less_gles[i,j] = gles[index,j]
                    less_gres[i,j] = gres[index,j]

            self.left.add_grids(nless, less_gles, less_gres,
                      l_ids, rank, size)

        if ngreater > 0:
            greater_gles = cvarray(format="d", shape=(ngreater,3), itemsize=sizeof(np.float64_t))
            greater_gres = cvarray(format="d", shape=(ngreater,3), itemsize=sizeof(np.float64_t))
            g_ids = cvarray(format="q", shape=(ngreater,), itemsize=sizeof(np.int64_t))

            for i in range(ngreater):
                index = greater_ids[i]
                g_ids[i] = gids[index]
                for j in range(3):
                    greater_gles[i,j] = gles[index,j]
                    greater_gres[i,j] = gres[index,j]

            self.right.add_grids(ngreater, greater_gles, greater_gres,
                      g_ids, rank, size)

        return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int should_i_split(self, int rank, int size):
        if self.node_id < size and self.node_id > 0:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void insert_grids(self,
                           int ngrids,
                           np.float64_t[:,:] gles,
                           np.float64_t[:,:] gres,
                           np.int64_t[:] gids,
                           int rank,
                           int size):

        if not should_i_build(self, rank, size) or ngrids == 0:
            return
        cdef int contained = 1
        cdef int check

        if ngrids == 1:
            # If we should continue to split based on parallelism, do so!
            if self.should_i_split(rank, size):
                self.geo_split(gles[0,:], gres[0,:], gids[0], rank, size)
                return

            for i in range(3):
                contained *= gles[0,i] <= self.left_edge[i]
                contained *= gres[0,i] >= self.right_edge[i]

            if contained == 1:
                # print 'Node fully contained, setting to grid: %i' % gids[0]
                self.grid = gids[0]
                assert(self.grid != -1)
                return

        # Split the grids
        check = self.split_grids(ngrids, gles, gres, gids, rank, size)
        # If check is -1, then we have found a place where there are no choices.
        # Exit out and set the node to None.
        if check == -1:
            self.grid = -1
        return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef split_grid(self,
                   np.float64_t[:] gle,
                   np.float64_t[:] gre,
                   int gid,
                   int rank,
                   int size):

        cdef int j
        cdef np.uint8_t[:] less_ids, greater_ids
        data = cvarray(format="d", shape=(1,2,3), itemsize=sizeof(np.float64_t))
        for j in range(3):
            data[0,0,j] = gle[j]
            data[0,1,j] = gre[j]

        less_ids = cvarray(format="B", shape=(1,), itemsize=sizeof(np.uint8_t))
        greater_ids = cvarray(format="B", shape=(1,), itemsize=sizeof(np.uint8_t))

        best_dim, split_pos, nless, ngreater = \
            kdtree_get_choices(1, data, self.left_edge, self.right_edge,
                              less_ids, greater_ids)

        # If best_dim is -1, then we have found a place where there are no choices.
        # Exit out and set the node to None.
        if best_dim == -1:
            return -1


        split = <Split *> malloc(sizeof(Split))
        split.dim = best_dim
        split.pos = split_pos

        # Create a Split
        self.divide(split)

        # Populate Left Node
        #print 'Inserting left node', self.left_edge, self.right_edge
        if nless == 1:
            self.left.insert_grid(gle, gre,
                         gid, rank, size)

        # Populate Right Node
        #print 'Inserting right node', self.left_edge, self.right_edge
        if ngreater == 1:
            self.right.insert_grid(gle, gre,
                         gid, rank, size)

        return 0

    #@cython.boundscheck(False)
    #@cython.wraparound(False)
    #@cython.cdivision(True)
    cdef int split_grids(self,
                           int ngrids,
                           np.float64_t[:,:] gles,
                           np.float64_t[:,:] gres,
                           np.int64_t[:] gids,
                           int rank,
                           int size):
        # Find a Split
        cdef int i, j, index
        cdef np.float64_t[:,:] less_gles, less_gres, greater_gles, greater_gres
        cdef np.int64_t[:] l_ids, g_ids
        if ngrids == 0: return 0

        data = cvarray(format="d", shape=(ngrids,2,3), itemsize=sizeof(np.float64_t))

        for i in range(ngrids):
            for j in range(3):
                data[i,0,j] = gles[i,j]
                data[i,1,j] = gres[i,j]

        less_ids = cvarray(format="B", shape=(ngrids,), itemsize=sizeof(np.uint8_t))
        greater_ids = cvarray(format="B", shape=(ngrids,), itemsize=sizeof(np.uint8_t))

        best_dim, split_pos, nless, ngreater = \
            kdtree_get_choices(ngrids, data, self.left_edge, self.right_edge,
                              less_ids, greater_ids)


        # If best_dim is -1, then we have found a place where there are no choices.
        # Exit out and set the node to None.
        if best_dim == -1:
            print 'Failed to split grids.'
            return -1

        split = <Split *> malloc(sizeof(Split))
        split.dim = best_dim
        split.pos = split_pos

        # Create a Split
        self.divide(split)

        less_index = cvarray(format="q", shape=(ngrids,), itemsize=sizeof(np.int64_t))
        greater_index = cvarray(format="q", shape=(ngrids,), itemsize=sizeof(np.int64_t))

        nless = 0
        ngreater = 0
        for i in range(ngrids):
            if less_ids[i] == 1:
                less_index[nless] = i
                nless += 1

            if greater_ids[i] == 1:
                greater_index[ngreater] = i
                ngreater += 1

        if nless > 0:
            less_gles = cvarray(format="d", shape=(nless,3), itemsize=sizeof(np.float64_t))
            less_gres = cvarray(format="d", shape=(nless,3), itemsize=sizeof(np.float64_t))
            l_ids = cvarray(format="q", shape=(nless,), itemsize=sizeof(np.int64_t))

            for i in range(nless):
                index = less_index[i]
                l_ids[i] = gids[index]
                for j in range(3):
                    less_gles[i,j] = gles[index,j]
                    less_gres[i,j] = gres[index,j]

            # Populate Left Node
            #print 'Inserting left node', self.left_edge, self.right_edge
            self.left.insert_grids(nless, less_gles, less_gres,
                         l_ids, rank, size)

        if ngreater > 0:
            greater_gles = cvarray(format="d", shape=(ngreater,3), itemsize=sizeof(np.float64_t))
            greater_gres = cvarray(format="d", shape=(ngreater,3), itemsize=sizeof(np.float64_t))
            g_ids = cvarray(format="q", shape=(ngreater,), itemsize=sizeof(np.int64_t))

            for i in range(ngreater):
                index = greater_index[i]
                g_ids[i] = gids[index]
                for j in range(3):
                    greater_gles[i,j] = gles[index,j]
                    greater_gres[i,j] = gres[index,j]

            # Populate Right Node
            #print 'Inserting right node', self.left_edge, self.right_edge
            self.right.insert_grids(ngreater, greater_gles, greater_gres,
                         g_ids, rank, size)

        return 0

    cdef geo_split(self,
                   np.float64_t[:] gle,
                   np.float64_t[:] gre,
                   int grid_id,
                   int rank,
                   int size):
        cdef int big_dim = 0
        cdef int i
        cdef np.float64_t v, my_max = 0.0

        for i in range(3):
            v = gre[i] - gle[i]
            if v > my_max:
                my_max = v
                big_dim = i

        new_pos = (gre[big_dim] + gle[big_dim])/2.

        lnew_gle = cvarray(format="d", shape=(3,), itemsize=sizeof(np.float64_t))
        lnew_gre = cvarray(format="d", shape=(3,), itemsize=sizeof(np.float64_t))
        rnew_gle = cvarray(format="d", shape=(3,), itemsize=sizeof(np.float64_t))
        rnew_gre = cvarray(format="d", shape=(3,), itemsize=sizeof(np.float64_t))

        for j in range(3):
            lnew_gle[j] = gle[j]
            lnew_gre[j] = gre[j]
            rnew_gle[j] = gle[j]
            rnew_gre[j] = gre[j]

        split = <Split *> malloc(sizeof(Split))
        split.dim = big_dim
        split.pos = new_pos

        # Create a Split
        self.divide(split)

        #lnew_gre[big_dim] = new_pos
        # Populate Left Node
        #print 'Inserting left node', self.left_edge, self.right_edge
        self.left.insert_grid(lnew_gle, lnew_gre,
                grid_id, rank, size)

        #rnew_gle[big_dim] = new_pos
        # Populate Right Node
        #print 'Inserting right node', self.left_edge, self.right_edge
        self.right.insert_grid(rnew_gle, rnew_gre,
                grid_id, rank, size)
        return

    cdef void divide(self, Split * split):
        # Create a Split
        self.split = split

        cdef np.float64_t[:] le = np.empty(3, dtype='float64')
        cdef np.float64_t[:] re = np.empty(3, dtype='float64')

        cdef int i
        for i in range(3):
            le[i] = self.left_edge[i]
            re[i] = self.right_edge[i]
        re[split.dim] = split.pos

        self.left = Node(self, None, None,
                         le, re, self.grid,
                         _lchild_id(self.node_id))

        re[split.dim] = self.right_edge[split.dim]
        le[split.dim] = split.pos
        self.right = Node(self, None, None,
                          le, re, self.grid,
                          _rchild_id(self.node_id))

        return
    #
    def kd_sum_volume(self):
        cdef np.float64_t vol = 1.0
        if (self.left is None) and (self.right is None):
            if self.grid == -1:
                return 0.0
            for i in range(3):
                vol *= self.right_edge[i] - self.left_edge[i]
            return vol
        else:
            return self.left.kd_sum_volume() + self.right.kd_sum_volume()

    def kd_node_check(self):
        assert (self.left is None) == (self.right is None)
        if (self.left is None) and (self.right is None):
            if self.grid != -1:
                return np.prod(self.right_edge - self.left_edge)
            else: return 0.0
        else:
            return self.left.kd_node_check()+self.right.kd_node_check()

    def kd_is_leaf(self):
        cdef int has_l_child = self.left == None
        cdef int has_r_child = self.right == None
        assert has_l_child == has_r_child
        return has_l_child

    cdef int _kd_is_leaf(self):
        if self.left is None or self.right is None:
            return 1
        return 0

    def depth_traverse(self, max_node=None):
        '''
        Yields a depth-first traversal of the kd tree always going to
        the left child before the right.
        '''
        current = self
        previous = None
        if max_node is None:
            max_node = np.inf
        while current is not None:
            yield current
            current, previous = step_depth(current, previous)
            if current is None: break
            if current.node_id >= max_node:
                current = current.parent
                previous = current.right

    def depth_first_touch(self, max_node=None):
        '''
        Yields a depth-first traversal of the kd tree always going to
        the left child before the right.
        '''
        current = self
        previous = None
        if max_node is None:
            max_node = np.inf
        while current is not None:
            if previous is None or previous.parent != current:
                yield current
            current, previous = step_depth(current, previous)
            if current is None: break
            if current.node_id >= max_node:
                current = current.parent
                previous = current.right

    def breadth_traverse(self):
        '''
        Yields a breadth-first traversal of the kd tree always going to
        the left child before the right.
        '''
        current = self
        previous = None
        while current is not None:
            yield current
            current, previous = step_depth(current, previous)


    def viewpoint_traverse(self, viewpoint):
        '''
        Yields a viewpoint dependent traversal of the kd-tree.  Starts
        with nodes furthest away from viewpoint.
        '''

        current = self
        previous = None
        while current is not None:
            yield current
            current, previous = step_viewpoint(current, previous, viewpoint)

    cdef int point_in_node(self,
                           np.float64_t[:] point):
        cdef int i
        cdef int inside = 1
        for i in range(3):
            inside *= self.left_edge[i] <= point[i]
            inside *= self.right_edge[i] > point[i]
        return inside

    cdef Node _find_node(self, np.float64_t[:] point):
        while self._kd_is_leaf() == 0:
            if point[self.split.dim] < self.split.pos:
                self = self.left
            else:
                self = self.right
        return self

    def find_node(self,
                  np.float64_t[:] point):
        """
        Find the AMRKDTree node enclosing a position
        """
        assert(self.point_in_node(point))
        return self._find_node(point)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline np.int64_t _lchild_id(np.int64_t node_id):
    return (node_id<<1)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline np.int64_t _rchild_id(np.int64_t node_id):
    return (node_id<<1) + 1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline np.int64_t _parent_id(np.int64_t node_id):
    return (node_id-1) >> 1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int should_i_build(Node node, int rank, int size):
    if (node.node_id < size) or (node.node_id >= 2*size):
        return 1
    elif node.node_id - size == rank:
        return 1
    else:
        return 0

def step_depth(Node current, Node previous):
    '''
    Takes a single step in the depth-first traversal
    '''
    if current._kd_is_leaf() == 1: # At a leaf, move back up
        previous = current
        current = current.parent

    elif current.parent is previous: # Moving down, go left first
        previous = current
        if current.left is not None:
            current = current.left
        elif current.right is not None:
            current = current.right
        else:
            current = current.parent

    elif current.left is previous: # Moving up from left, go right
        previous = current
        if current.right is not None:
            current = current.right
        else:
            current = current.parent

    elif current.right is previous: # Moving up from right child, move up
        previous = current
        current = current.parent

    return current, previous

def step_viewpoint(Node current,
                   Node previous,
                   viewpoint):
    '''
    Takes a single step in the viewpoint based traversal.  Always
    goes to the node furthest away from viewpoint first.
    '''
    if current._kd_is_leaf() == 1: # At a leaf, move back up
        previous = current
        current = current.parent
    elif current.split.dim is None: # This is a dead node
        previous = current
        current = current.parent

    elif current.parent is previous: # Moving down
        previous = current
        if viewpoint[current.split.dim] <= current.split.pos:
            if current.right is not None:
                current = current.right
            else:
                previous = current.right
        else:
            if current.left is not None:
                current = current.left
            else:
                previous = current.left

    elif current.right is previous: # Moving up from right
        previous = current
        if viewpoint[current.split.dim] <= current.split.pos:
            if current.left is not None:
                current = current.left
            else:
                current = current.parent
        else:
            current = current.parent

    elif current.left is previous: # Moving up from left child
        previous = current
        if viewpoint[current.split.dim] > current.split.pos:
            if current.right is not None:
                current = current.right
            else:
                current = current.parent
        else:
            current = current.parent

    return current, previous

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef kdtree_get_choices(int n_grids,
                        np.float64_t[:,:,:] data,
                        np.float64_t[:] l_corner,
                        np.float64_t[:] r_corner,
                        np.uint8_t[:] less_ids,
                        np.uint8_t[:] greater_ids,
                       ):
    cdef int i, j, k, dim, n_unique, best_dim, my_split
    cdef np.float64_t split
    cdef np.float64_t[:,:] uniquedims
    cdef np.float64_t[:] uniques
    uniquedims = cvarray(format="d", shape=(3, 2*n_grids), itemsize=sizeof(np.float64_t))
    my_max = 0
    my_split = 0
    best_dim = -1
    for dim in range(3):
        n_unique = 0
        uniques = uniquedims[dim]
        for i in range(n_grids):
            # Check for disqualification
            for j in range(2):
                # print "Checking against", i,j,dim,data[i,j,dim]
                if not (l_corner[dim] < data[i][j][dim] and
                        data[i][j][dim] < r_corner[dim]):
                    # print "Skipping ", data[i,j,dim], l_corner[dim], r_corner[dim]
                    continue
                skipit = 0
                # Add our left ...
                for k in range(n_unique):
                    if uniques[k] == data[i][j][dim]:
                        skipit = 1
                        # print "Identified", uniques[k], data[i,j,dim], n_unique
                        break
                if skipit == 0:
                    uniques[n_unique] = data[i][j][dim]
                    n_unique += 1
        if n_unique > my_max:
            best_dim = dim
            my_max = n_unique
            my_split = (n_unique-1)/2
    if best_dim == -1:
        return -1, 0, 0, 0
    # I recognize how lame this is.
    cdef np.ndarray[np.float64_t, ndim=1] tarr = np.empty(my_max, dtype='float64')
    for i in range(my_max):
        # print "Setting tarr: ", i, uniquedims[best_dim][i]
        tarr[i] = uniquedims[best_dim][i]
    tarr.sort()
    split = tarr[my_split]
    cdef int nless=0, ngreater=0
    for i in range(n_grids):
        if data[i][0][best_dim] < split:
            less_ids[i] = 1
            nless += 1
        else:
            less_ids[i] = 0
        if data[i][1][best_dim] > split:
            greater_ids[i] = 1
            ngreater += 1
        else:
            greater_ids[i] = 0

    # Return out unique values
    return best_dim, split, nless, ngreater
