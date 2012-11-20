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
    cdef np.float64_t * xp
    cdef np.float64_t * yp
    cdef np.float64_t * zp
    cdef GridTree tree
    cdef np.int64_t * point_grids

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

        self.xp = <np.float64_t *> malloc(sizeof(np.float64_t) *
                                          num_points)
        self.yp = <np.float64_t *> malloc(sizeof(np.float64_t) *
                                          num_points)
        self.zp = <np.float64_t *> malloc(sizeof(np.float64_t) *
                                          num_points)
        self.point_grids = <np.int64_t *> malloc(sizeof(np.int64_t) *
                                              num_points)
        
        for i in range(num_points) :
            self.xp[i] = x[i]
            self.yp[i] = y[i]
            self.zp[i] = z[i]
            self.point_grids[i] = -1
            
        self.tree = tree

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def find_points_in_tree(self) :

        cdef np.ndarray[np.int64_t, ndim=1] pt_grids
        cdef int i, j

        pt_grids = np.zeros(self.num_points, dtype='int64')

        for i in range(self.num_points) :

            for j in range(self.tree.num_root_grids) :

                self.check_position(i, self.xp[i], self.yp[i], self.zp[i],
				    &self.tree.root_grids[j])

        for i in range(self.num_points) :
            pt_grids[i] = self.point_grids[i]
        
        return pt_grids

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void check_position(self,
			     np.int64_t pt_index, 
			     np.float64_t x,
			     np.float64_t y,
			     np.float64_t z,
			     GridTreeNode * grid) :

        cdef int i
        cdef np.uint8_t in_grid
	
        in_grid = self.is_in_grid(x, y, z, grid)

        if in_grid :

            if grid.num_children > 0 :

                for i in range(grid.num_children) :

                    self.check_position(pt_index, x, y, z, grid.children[i])

            else :

                self.point_grids[pt_index] = grid.index

        else :

            return
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.uint8_t is_in_grid(self,
			 np.float64_t x,
			 np.float64_t y,
			 np.float64_t z,
			 GridTreeNode * grid) :

        cdef np.uint8_t xcond, ycond, zcond, cond
            
        xcond = x >= grid.left_edge[0] and x < grid.right_edge[0]
        ycond = y >= grid.left_edge[1] and y < grid.right_edge[1]
        zcond = z >= grid.left_edge[2] and z < grid.right_edge[2]
	
        cond = xcond and ycond
        cond = cond and zcond

        return cond
