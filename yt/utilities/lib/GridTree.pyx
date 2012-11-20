"""
Matching points on the grid to specific grids

Author: John ZuHone <jzuhone@gmail.com>
Affiliation: NASA/Goddard Space Flight Center
Homepage: http://yt-project.org/
License:
Copyright (C) 2012 John ZuHone.  All Rights Reserved.

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
cimport cython

from libc.stdlib cimport malloc, free

cdef struct GridTreeNode :
    int num_children
    int level
    int index
    np.float64_t left_edge[3]
    np.float64_t right_edge[3]
    GridTreeNode **children
                
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef GridTreeNode Grid_initialize(np.ndarray[np.float64_t, ndim=1] le,
                                  np.ndarray[np.float64_t, ndim=1] re,
                                  int num_children, int level, int index) :

    cdef GridTreeNode node
    cdef int i

    node.index = index
    node.level = level
    for i in range(3) :
        node.left_edge[i] = le[i]
        node.right_edge[i] = re[i]
    node.num_children = num_children
    
    if num_children > 0:
        node.children = <GridTreeNode **> malloc(sizeof(GridTreeNode *) *
                                                 num_children)
        for i in range(num_children) :
            node.children[i] = NULL
    else :
        node.children = NULL

    return node

cdef class GridTree :

    cdef GridTreeNode *grids
    cdef GridTreeNode *root_grids
    cdef int num_grids
    cdef int num_root_grids
    cdef int num_leaf_grids
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __cinit__(self, int num_grids, 
                  np.ndarray[np.float64_t, ndim=2] left_edge,
                  np.ndarray[np.float64_t, ndim=2] right_edge,
                  np.ndarray[np.int64_t, ndim=1] parent_ind,
                  np.ndarray[np.int64_t, ndim=1] level,
                  np.ndarray[np.int64_t, ndim=1] num_children) :

        cdef int i, j, k
        cdef np.ndarray[np.int64_t, ndim=1] child_ptr

        child_ptr = np.zeros(num_grids, dtype='int64')

        self.num_grids = num_grids
        self.num_root_grids = 0
        self.num_leaf_grids = 0
        
        self.grids = <GridTreeNode *> malloc(sizeof(GridTreeNode) *
                                             num_grids)
                
        for i in range(num_grids) :

            self.grids[i] = Grid_initialize(left_edge[i,:],
                                            right_edge[i,:],
                                            num_children[i],
                                            level[i], i)
            if level[i] == 0 :
                self.num_root_grids += 1

            if num_children[i] == 0 : self.num_leaf_grids += 1

        self.root_grids = <GridTreeNode *> malloc(sizeof(GridTreeNode) *
                                                  self.num_root_grids)
                
        k = 0
        
        for i in range(num_grids) :

            j = parent_ind[i]
            
            if j >= 0:
                
                self.grids[j].children[child_ptr[j]] = &self.grids[i]

                child_ptr[j] += 1

            else :

                self.root_grids[k] = self.grids[i] 
                
                k = k + 1
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def return_tree_info(self) :

        cdef int i, j
        
        levels = []
        indices = []
        nchild = []
        children = []
        
        for i in range(self.num_grids) : 

            childs = []
            
            levels.append(self.grids[i].level)
            indices.append(self.grids[i].index)
            nchild.append(self.grids[i].num_children)
            for j in range(self.grids[i].num_children) :
                childs.append(self.grids[i].children[j].index)
            children.append(childs)

        return indices, levels, nchild, children
    
cdef class MatchPointsToGrids :

    cdef int num_points
    cdef int num_idxs_left
    cdef int num_grids_walked
    cdef np.float64_t * xp
    cdef np.float64_t * yp
    cdef np.float64_t * zp
    cdef GridTree tree
    cdef np.int64_t * all_idxs
    cdef np.int64_t * point_grids
    cdef np.int64_t * in_grid

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __cinit__(self, GridTree tree,
                  int num_points, 
                  np.ndarray[np.float64_t, ndim=1] x,
                  np.ndarray[np.float64_t, ndim=1] y,
                  np.ndarray[np.float64_t, ndim=1] z) :

        cdef int i
        
        self.num_points = num_points
        self.num_idxs_left = num_points

        self.xp = <np.float64_t *> malloc(sizeof(np.float64_t) *
                                          num_points)
        self.yp = <np.float64_t *> malloc(sizeof(np.float64_t) *
                                          num_points)
        self.zp = <np.float64_t *> malloc(sizeof(np.float64_t) *
                                          num_points)
        self.all_idxs = <np.int64_t *> malloc(sizeof(np.int64_t) *
                                              num_points)
        self.in_grid = <np.int64_t *> malloc(sizeof(np.int64_t) *
                                             num_points)
        self.point_grids = <np.int64_t *> malloc(sizeof(np.int64_t) *
                                              num_points)
        
        for i in range(num_points) :
            self.xp[i] = x[i]
            self.yp[i] = y[i]
            self.zp[i] = z[i]
            self.all_idxs[i] = i
            self.in_grid[i] = 0
            self.point_grids[i] = -1
            
        self.tree = tree

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def find_points_in_tree(self) :

        cdef np.ndarray[np.int64_t, ndim=1] pt_grids
        cdef int i

        pt_grids = np.zeros(self.num_points, dtype='int64')
        
        for i in range(self.tree.num_root_grids) :

            self.check_positions(&self.tree.root_grids[i])

        for i in range(self.num_points) :
            pt_grids[i] = self.point_grids[i]
        
        return pt_grids

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void check_positions(self, GridTreeNode * grid) :

        cdef int i
        
        if self.num_idxs_left == 0 : return

        if grid.num_children > 0 :

            for i in range(grid.num_children) :
                
                self.check_positions(grid.children[i])
                
        else :

            self.is_in_grid(grid)

            for i in range(self.num_points) :

                if self.in_grid[i] == 1 :

                    #print self.num_idxs_left
                    
                    self.point_grids[i] = grid.index
                    self.all_idxs[i] = -1
                    self.num_idxs_left -= 1
                    
            return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void is_in_grid(self, GridTreeNode * grid) :

        cdef int i
        cdef np.uint8_t xcond, ycond, zcond, cond
        
        for i in range(self.num_points) :

            if self.all_idxs[i] >= 0 :
            
                xcond = self.xp[i] >= grid.left_edge[0] and self.xp[i] < grid.right_edge[0]
                ycond = self.yp[i] >= grid.left_edge[1] and self.yp[i] < grid.right_edge[1]
                zcond = self.zp[i] >= grid.left_edge[2] and self.zp[i] < grid.right_edge[2]

                cond = xcond and ycond
                cond = cond and zcond
                
                if cond :
                    self.in_grid[i] = 1
                else :
                    self.in_grid[i] = 0

            else :

                self.in_grid[i] = 0
