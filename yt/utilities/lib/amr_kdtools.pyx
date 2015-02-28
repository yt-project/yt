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

cdef extern from "stdlib.h":
    # NOTE that size_t might not be int
    void *alloca(int)


DEF Nch = 4

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef class Node:

    def __cinit__(self,
                  Node parent,
                  Node left,
                  Node right,
                  np.ndarray[np.float64_t, ndim=1] left_edge,
                  np.ndarray[np.float64_t, ndim=1] right_edge,
                  int grid,
                  np.int64_t node_id):
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

    def get_left_edge(self):
        return get_left_edge(self)

    def get_right_edge(self):
        return get_right_edge(self)

    def set_left_edge(self, np.ndarray[np.float64_t, ndim=1] left_edge):
        cdef int i
        for i in range(3):
            self.left_edge[i] = left_edge[i]

    def set_right_edge(self, np.ndarray[np.float64_t, ndim=1] right_edge):
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

def get_left_edge(Node node):
    le = np.empty(3, dtype='float64')
    for i in range(3):
        le[i] = node.left_edge[i]
    return le

def get_right_edge(Node node):
    re = np.empty(3, dtype='float64')
    for i in range(3):
        re[i] = node.right_edge[i]
    return re

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

def kd_traverse(Node trunk, viewpoint=None):
    if viewpoint is None:
        for node in depth_traverse(trunk):
            if _kd_is_leaf(node) == 1 and node.grid != -1:
                yield node
    else:
        for node in viewpoint_traverse(trunk, viewpoint):
            if _kd_is_leaf(node) == 1 and node.grid != -1:
                yield node

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef add_grid(Node node,
                   np.float64_t *gle,
                   np.float64_t *gre,
                   int gid,
                   int rank,
                   int size):

    if not should_i_build(node, rank, size):
        return

    if _kd_is_leaf(node) == 1:
        insert_grid(node, gle, gre, gid, rank, size)
    else:
        less_id = gle[node.split.dim] < node.split.pos
        if less_id:
            add_grid(node.left, gle, gre,
                     gid, rank, size)

        greater_id = gre[node.split.dim] > node.split.pos
        if greater_id:
            add_grid(node.right, gle, gre,
                     gid, rank, size)
    return

def add_pygrid(Node node,
                   np.ndarray[np.float64_t, ndim=1] gle,
                   np.ndarray[np.float64_t, ndim=1] gre,
                   int gid,
                   int rank,
                   int size):

    """
    The entire purpose of this function is to move everything from ndarrays
    to internal C pointers.
    """
    pgles = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
    pgres = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
    cdef int j
    for j in range(3):
        pgles[j] = gle[j]
        pgres[j] = gre[j]

    add_grid(node, pgles, pgres, gid, rank, size)
    free(pgles)
    free(pgres)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef insert_grid(Node node,
                np.float64_t *gle,
                np.float64_t *gre,
                int grid_id,
                int rank,
                int size):
    if not should_i_build(node, rank, size):
        return

    # If we should continue to split based on parallelism, do so!
    if should_i_split(node, rank, size):
        geo_split(node, gle, gre, grid_id, rank, size)
        return

    cdef int contained = 1
    for i in range(3):
        if gle[i] > node.left_edge[i] or\
           gre[i] < node.right_edge[i]:
            contained *= 0

    if contained == 1:
        node.grid = grid_id
        assert(node.grid != -1)
        return

    # Split the grid
    cdef int check = split_grid(node, gle, gre, grid_id, rank, size)
    # If check is -1, then we have found a place where there are no choices.
    # Exit out and set the node to None.
    if check == -1:
        node.grid = -1
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def add_pygrids(Node node,
                    int ngrids,
                    np.ndarray[np.float64_t, ndim=2] gles,
                    np.ndarray[np.float64_t, ndim=2] gres,
                    np.ndarray[np.int64_t, ndim=1] gids,
                    int rank,
                    int size):
    """
    The entire purpose of this function is to move everything from ndarrays
    to internal C pointers.
    """
    pgles = <np.float64_t **> malloc(ngrids * sizeof(np.float64_t*))
    pgres = <np.float64_t **> malloc(ngrids * sizeof(np.float64_t*))
    pgids = <np.int64_t *> malloc(ngrids * sizeof(np.int64_t))
    for i in range(ngrids):
        pgles[i] = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
        pgres[i] = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
        pgids[i] = gids[i]
        for j in range(3):
            pgles[i][j] = gles[i, j]
            pgres[i][j] = gres[i, j]

    add_grids(node, ngrids, pgles, pgres, pgids, rank, size)

    for i in range(ngrids):
        free(pgles[i])
        free(pgres[i])
    free(pgles)
    free(pgres)
    free(pgids)



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef add_grids(Node node,
                    int ngrids,
                    np.float64_t **gles,
                    np.float64_t **gres,
                    np.int64_t *gids,
                    int rank,
                    int size):
    cdef int i, j, nless, ngreater
    cdef np.int64_t gid
    if not should_i_build(node, rank, size):
        return

    if _kd_is_leaf(node) == 1:
        insert_grids(node, ngrids, gles, gres, gids, rank, size)
        return

    less_ids= <np.int64_t *> malloc(ngrids * sizeof(np.int64_t))
    greater_ids = <np.int64_t *> malloc(ngrids * sizeof(np.int64_t))

    nless = 0
    ngreater = 0
    for i in range(ngrids):
        if gles[i][node.split.dim] < node.split.pos:
            less_ids[nless] = i
            nless += 1

        if gres[i][node.split.dim] > node.split.pos:
            greater_ids[ngreater] = i
            ngreater += 1

    #print 'nless: %i' % nless
    #print 'ngreater: %i' % ngreater

    less_gles = <np.float64_t **> malloc(nless * sizeof(np.float64_t*))
    less_gres = <np.float64_t **> malloc(nless * sizeof(np.float64_t*))
    l_ids = <np.int64_t *> malloc(nless * sizeof(np.int64_t))
    for i in range(nless):
        less_gles[i] = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
        less_gres[i] = <np.float64_t *> malloc(3 * sizeof(np.float64_t))

    greater_gles = <np.float64_t **> malloc(ngreater * sizeof(np.float64_t*))
    greater_gres = <np.float64_t **> malloc(ngreater * sizeof(np.float64_t*))
    g_ids = <np.int64_t *> malloc(ngreater * sizeof(np.int64_t))
    for i in range(ngreater):
        greater_gles[i] = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
        greater_gres[i] = <np.float64_t *> malloc(3 * sizeof(np.float64_t))

    cdef int index
    for i in range(nless):
        index = less_ids[i]
        l_ids[i] = gids[index]
        for j in range(3):
            less_gles[i][j] = gles[index][j]
            less_gres[i][j] = gres[index][j]

    if nless > 0:
        add_grids(node.left, nless, less_gles, less_gres,
                  l_ids, rank, size)

    for i in range(ngreater):
        index = greater_ids[i]
        g_ids[i] = gids[index]
        for j in range(3):
            greater_gles[i][j] = gles[index][j]
            greater_gres[i][j] = gres[index][j]

    if ngreater > 0:
        add_grids(node.right, ngreater, greater_gles, greater_gres,
                  g_ids, rank, size)

    for i in range(nless):
        free(less_gles[i])
        free(less_gres[i])
    free(less_gles)
    free(less_gres)
    free(less_ids)
    free(l_ids)
    for i in range(ngreater):
        free(greater_gles[i])
        free(greater_gres[i])
    free(greater_gles)
    free(greater_gres)
    free(greater_ids)
    free(g_ids)

    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int should_i_split(Node node, int rank, int size):
    if node.node_id < size and node.node_id > 0:
        return 1
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void insert_grids(Node node,
                       int ngrids,
                       np.float64_t **gles,
                       np.float64_t **gres,
                       np.int64_t *gids,
                       int rank,
                       int size):

    if not should_i_build(node, rank, size) or ngrids == 0:
        return
    cdef int contained = 1
    cdef int check

    if ngrids == 1:
        # If we should continue to split based on parallelism, do so!
        if should_i_split(node, rank, size):
            geo_split(node, gles[0], gres[0], gids[0], rank, size)
            return

        for i in range(3):
            contained *= gles[0][i] <= node.left_edge[i]
            contained *= gres[0][i] >= node.right_edge[i]

        if contained == 1:
            # print 'Node fully contained, setting to grid: %i' % gids[0]
            node.grid = gids[0]
            assert(node.grid != -1)
            return

    # Split the grids
    check = split_grids(node, ngrids, gles, gres, gids, rank, size)
    # If check is -1, then we have found a place where there are no choices.
    # Exit out and set the node to None.
    if check == -1:
        node.grid = -1
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef split_grid(Node node,
               np.float64_t *gle,
               np.float64_t *gre,
               int gid,
               int rank,
               int size):

    cdef int j
    data = <np.float64_t ***> alloca(sizeof(np.float64_t**))
    data[0] = <np.float64_t **> alloca(2 * sizeof(np.float64_t*))
    for j in range(2):
        data[0][j] = <np.float64_t *> alloca(3 * sizeof(np.float64_t))
    for j in range(3):
        data[0][0][j] = gle[j]
        data[0][1][j] = gre[j]

    less_ids = <np.uint8_t *> alloca(1 * sizeof(np.uint8_t))
    greater_ids = <np.uint8_t *> alloca(1 * sizeof(np.uint8_t))

    best_dim, split_pos, nless, ngreater = \
        kdtree_get_choices(1, data, node.left_edge, node.right_edge,
                          less_ids, greater_ids)

    # If best_dim is -1, then we have found a place where there are no choices.
    # Exit out and set the node to None.
    if best_dim == -1:
        print 'Failed to split grid.'
        return -1


    split = <Split *> malloc(sizeof(Split))
    split.dim = best_dim
    split.pos = split_pos

    # Create a Split
    divide(node, split)

    # Populate Left Node
    #print 'Inserting left node', node.left_edge, node.right_edge
    if nless == 1:
        insert_grid(node.left, gle, gre,
                     gid, rank, size)

    # Populate Right Node
    #print 'Inserting right node', node.left_edge, node.right_edge
    if ngreater == 1:
        insert_grid(node.right, gle, gre,
                     gid, rank, size)

    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef kdtree_get_choices(int n_grids,
                        np.float64_t ***data,
                        np.float64_t *l_corner,
                        np.float64_t *r_corner,
                        np.uint8_t *less_ids,
                        np.uint8_t *greater_ids,
                       ):
    cdef int i, j, k, dim, n_unique, best_dim, n_best, addit, my_split
    cdef np.float64_t split
    cdef np.float64_t **uniquedims
    cdef np.float64_t *uniques
    uniquedims = <np.float64_t **> malloc(3 * sizeof(np.float64_t*))
    for i in range(3):
        uniquedims[i] = <np.float64_t *> \
                malloc(2*n_grids * sizeof(np.float64_t))
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

    for i in range(3):
        free(uniquedims[i])
    free(uniquedims)

    # Return out unique values
    return best_dim, split, nless, ngreater

#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.cdivision(True)
cdef int split_grids(Node node,
                       int ngrids,
                       np.float64_t **gles,
                       np.float64_t **gres,
                       np.int64_t *gids,
                       int rank,
                       int size):
    # Find a Split
    cdef int i, j, k

    data = <np.float64_t ***> malloc(ngrids * sizeof(np.float64_t**))
    for i in range(ngrids):
        data[i] = <np.float64_t **> malloc(2 * sizeof(np.float64_t*))
        for j in range(2):
            data[i][j] = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
        for j in range(3):
            data[i][0][j] = gles[i][j]
            data[i][1][j] = gres[i][j]

    less_ids = <np.uint8_t *> malloc(ngrids * sizeof(np.uint8_t))
    greater_ids = <np.uint8_t *> malloc(ngrids * sizeof(np.uint8_t))

    best_dim, split_pos, nless, ngreater = \
        kdtree_get_choices(ngrids, data, node.left_edge, node.right_edge,
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
    divide(node, split)

    less_index = <np.int64_t *> malloc(ngrids * sizeof(np.int64_t))
    greater_index = <np.int64_t *> malloc(ngrids * sizeof(np.int64_t))

    nless = 0
    ngreater = 0
    for i in range(ngrids):
        if less_ids[i] == 1:
            less_index[nless] = i
            nless += 1

        if greater_ids[i] == 1:
            greater_index[ngreater] = i
            ngreater += 1

    less_gles = <np.float64_t **> malloc(nless * sizeof(np.float64_t*))
    less_gres = <np.float64_t **> malloc(nless * sizeof(np.float64_t*))
    l_ids = <np.int64_t *> malloc(nless * sizeof(np.int64_t))
    for i in range(nless):
        less_gles[i] = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
        less_gres[i] = <np.float64_t *> malloc(3 * sizeof(np.float64_t))

    greater_gles = <np.float64_t **> malloc(ngreater * sizeof(np.float64_t*))
    greater_gres = <np.float64_t **> malloc(ngreater * sizeof(np.float64_t*))
    g_ids = <np.int64_t *> malloc(ngreater * sizeof(np.int64_t))
    for i in range(ngreater):
        greater_gles[i] = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
        greater_gres[i] = <np.float64_t *> malloc(3 * sizeof(np.float64_t))

    cdef int index
    for i in range(nless):
        index = less_index[i]
        l_ids[i] = gids[index]
        for j in range(3):
            less_gles[i][j] = gles[index][j]
            less_gres[i][j] = gres[index][j]

    if nless > 0:
        # Populate Left Node
        #print 'Inserting left node', node.left_edge, node.right_edge
        insert_grids(node.left, nless, less_gles, less_gres,
                     l_ids, rank, size)

    for i in range(ngreater):
        index = greater_index[i]
        g_ids[i] = gids[index]
        for j in range(3):
            greater_gles[i][j] = gles[index][j]
            greater_gres[i][j] = gres[index][j]

    if ngreater > 0:
        # Populate Right Node
        #print 'Inserting right node', node.left_edge, node.right_edge
        insert_grids(node.right, ngreater, greater_gles, greater_gres,
                     g_ids, rank, size)

    for i in range(nless):
        free(less_gles[i])
        free(less_gres[i])
    free(less_gles)
    free(less_gres)
    free(less_ids)
    free(less_index)
    free(l_ids)
    for i in range(ngreater):
        free(greater_gles[i])
        free(greater_gres[i])
    free(greater_gles)
    free(greater_gres)
    free(greater_ids)
    free(greater_index)
    free(g_ids)

    for i in range(ngrids):
        for j in range(2):
            free(data[i][j])
        free(data[i])
    free(data)

    return 0

cdef geo_split(Node node,
               np.float64_t *gle,
               np.float64_t *gre,
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

    lnew_gle = <np.float64_t *> alloca(3 * sizeof(np.float64_t))
    lnew_gre = <np.float64_t *> alloca(3 * sizeof(np.float64_t))
    rnew_gle = <np.float64_t *> alloca(3 * sizeof(np.float64_t))
    rnew_gre = <np.float64_t *> alloca(3 * sizeof(np.float64_t))

    for j in range(3):
        lnew_gle[j] = gle[j]
        lnew_gre[j] = gre[j]
        rnew_gle[j] = gle[j]
        rnew_gre[j] = gre[j]

    split = <Split *> malloc(sizeof(Split))
    split.dim = big_dim
    split.pos = new_pos

    # Create a Split
    divide(node, split)

    #lnew_gre[big_dim] = new_pos
    # Populate Left Node
    #print 'Inserting left node', node.left_edge, node.right_edge
    insert_grid(node.left, lnew_gle, lnew_gre,
            grid_id, rank, size)

    #rnew_gle[big_dim] = new_pos
    # Populate Right Node
    #print 'Inserting right node', node.left_edge, node.right_edge
    insert_grid(node.right, rnew_gle, rnew_gre,
            grid_id, rank, size)
    return

cdef void divide(Node node, Split * split):
    # Create a Split
    node.split = split

    cdef np.ndarray[np.float64_t, ndim=1] le = np.zeros(3, dtype='float64')
    cdef np.ndarray[np.float64_t, ndim=1] re = np.zeros(3, dtype='float64')

    cdef int i
    for i in range(3):
        le[i] = node.left_edge[i]
        re[i] = node.right_edge[i]
    re[split.dim] = split.pos

    node.left = Node(node, None, None,
                     le, re, node.grid,
                     _lchild_id(node.node_id))

    re[split.dim] = node.right_edge[split.dim]
    le[split.dim] = split.pos
    node.right = Node(node, None, None,
                      le, re, node.grid,
                      _rchild_id(node.node_id))

    return
#
def kd_sum_volume(Node node):
    cdef np.float64_t vol = 1.0
    if (node.left is None) and (node.right is None):
        if node.grid == -1:
            return 0.0
        for i in range(3):
            vol *= node.right_edge[i] - node.left_edge[i]
        return vol
    else:
        return kd_sum_volume(node.left) + kd_sum_volume(node.right)

def kd_node_check(Node node):
    assert (node.left is None) == (node.right is None)
    if (node.left is None) and (node.right is None):
        if node.grid != -1:
            return np.prod(node.right_edge - node.left_edge)
        else: return 0.0
    else:
        return kd_node_check(node.left)+kd_node_check(node.right)

def kd_is_leaf(Node node):
    cdef int has_l_child = node.left == None
    cdef int has_r_child = node.right == None
    assert has_l_child == has_r_child
    return has_l_child

cdef int _kd_is_leaf(Node node):
    if node.left is None or node.right is None:
        return 1
    return 0

def step_depth(Node current, Node previous):
    '''
    Takes a single step in the depth-first traversal
    '''
    if _kd_is_leaf(current) == 1: # At a leaf, move back up
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

def depth_traverse(Node trunk, max_node=None):
    '''
    Yields a depth-first traversal of the kd tree always going to
    the left child before the right.
    '''
    current = trunk
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

def depth_first_touch(Node tree, max_node=None):
    '''
    Yields a depth-first traversal of the kd tree always going to
    the left child before the right.
    '''
    current = tree
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

def breadth_traverse(Node tree):
    '''
    Yields a breadth-first traversal of the kd tree always going to
    the left child before the right.
    '''
    current = tree
    previous = None
    while current is not None:
        yield current
        current, previous = step_depth(current, previous)


def viewpoint_traverse(Node tree, viewpoint):
    '''
    Yields a viewpoint dependent traversal of the kd-tree.  Starts
    with nodes furthest away from viewpoint.
    '''

    current = tree
    previous = None
    while current is not None:
        yield current
        current, previous = step_viewpoint(current, previous, viewpoint)

def step_viewpoint(Node current,
                   Node previous,
                   viewpoint):
    '''
    Takes a single step in the viewpoint based traversal.  Always
    goes to the node furthest away from viewpoint first.
    '''
    if _kd_is_leaf(current) == 1: # At a leaf, move back up
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

cdef int point_in_node(Node node,
                       np.ndarray[np.float64_t, ndim=1] point):
    cdef int i
    cdef int inside = 1
    for i in range(3):
        inside *= node.left_edge[i] <= point[i]
        inside *= node.right_edge[i] > point[i]
    return inside

cdef Node _find_node(Node node, np.float64_t *point):
    while _kd_is_leaf(node) == 0:
        if point[node.split.dim] < node.split.pos:
            node = node.left
        else:
            node = node.right
    return node

def find_node(Node node,
              np.ndarray[np.float64_t, ndim=1] point):
    """
    Find the AMRKDTree node enclosing a position
    """
    assert(point_in_node(node, point))
    return _find_node(node, <np.float64_t *> point.data)

