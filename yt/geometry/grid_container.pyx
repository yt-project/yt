"""
Matching points on the grid to specific grids



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
from libc.math cimport rint
from yt.utilities.lib.bitarray cimport bitarray

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef GridTreeNode Grid_initialize(np.ndarray[np.float64_t, ndim=1] le,
                                  np.ndarray[np.float64_t, ndim=1] re,
                                  np.ndarray[np.int32_t, ndim=1] dims,
                                  int num_children, int level, int index):

    cdef GridTreeNode node
    cdef int i

    node.index = index
    node.level = level
    for i in range(3):
        node.left_edge[i] = le[i]
        node.right_edge[i] = re[i]
        node.dims[i] = dims[i]
        node.dds[i] = (re[i] - le[i])/dims[i]
        node.start_index[i] = <np.int64_t> rint(le[i] / node.dds[i])
    node.num_children = num_children
    if num_children <= 0:
        node.children = NULL
        return node
    node.children = <GridTreeNode **> malloc(
            sizeof(GridTreeNode *) * num_children)
    for i in range(num_children):
        node.children[i] = NULL

    return node

cdef class GridTree:
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __cinit__(self, int num_grids, 
                  np.ndarray[np.float64_t, ndim=2] left_edge,
                  np.ndarray[np.float64_t, ndim=2] right_edge,
                  np.ndarray[np.int32_t, ndim=2] dimensions,
                  np.ndarray[np.int64_t, ndim=1] parent_ind,
                  np.ndarray[np.int64_t, ndim=1] level,
                  np.ndarray[np.int64_t, ndim=1] num_children):

        cdef int i, j, k
        cdef np.ndarray[np.int64_t, ndim=1] child_ptr

        child_ptr = np.zeros(num_grids, dtype='int64')

        self.num_grids = num_grids
        self.num_root_grids = 0
        self.num_leaf_grids = 0
        
        self.grids = <GridTreeNode *> malloc(
                sizeof(GridTreeNode) * num_grids)
                
        for i in range(num_grids):
            self.grids[i] = Grid_initialize(left_edge[i,:],
                                            right_edge[i,:],
                                            dimensions[i,:],
                                            num_children[i],
                                            level[i], i)
            if level[i] == 0:
                self.num_root_grids += 1
            if num_children[i] == 0:
                self.num_leaf_grids += 1

        self.root_grids = <GridTreeNode *> malloc(
                sizeof(GridTreeNode) * self.num_root_grids)
        k = 0
        for i in range(num_grids):
            j = parent_ind[i]
            if j >= 0:
                self.grids[j].children[child_ptr[j]] = &self.grids[i]
                child_ptr[j] += 1
            else:
                if k >= self.num_root_grids:
                    raise RuntimeError
                self.root_grids[k] = self.grids[i] 
                k = k + 1

    def __iter__(self):
        yield self
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def return_tree_info(self):
        cdef int i, j
        levels = []
        indices = []
        nchild = []
        children = []
        for i in range(self.num_grids): 
            childs = []
            levels.append(self.grids[i].level)
            indices.append(self.grids[i].index)
            nchild.append(self.grids[i].num_children)
            for j in range(self.grids[i].num_children):
                childs.append(self.grids[i].children[j].index)
            children.append(childs)
        return indices, levels, nchild, children

    cdef void setup_data(self, GridVisitorData *data):
        data.index = 0
        data.global_index = 0
        data.n_tuples = 0
        data.child_tuples = NULL
        data.array = NULL
        data.ref_factor = 2 #### FIX THIS

    cdef void visit_grids(self, GridVisitorData *data,
                          grid_visitor_function *func,
                          SelectorObject selector):
        cdef int i, n
        # Because of confusion about mapping of children to parents, we are
        # going to do this the stupid way for now.
        cdef GridTreeNode *grid
        for i in range(self.num_root_grids):
            grid = &self.root_grids[i]
            self.recursively_visit_grid(data, func, selector, grid)
        grid_visitors.free_tuples(data)

    cdef void recursively_visit_grid(self, GridVisitorData *data,
                                     grid_visitor_function *func,
                                     SelectorObject selector,
                                     GridTreeNode *grid):
        cdef int i
        data.grid = grid
        if selector.select_bbox(grid.left_edge, grid.right_edge) == 0:
            # Note that this does not increment the global_index.
            return
        grid_visitors.setup_tuples(data)
        selector.visit_grid_cells(data, func)
        for i in range(grid.num_children):
            self.recursively_visit_grid(data, func, selector, grid.children[i])

    def count(self, SelectorObject selector):
        cdef GridVisitorData data
        self.setup_data(&data)
        cdef np.uint64_t size = 0
        data.array = <void*>(&size)
        self.visit_grids(&data,  grid_visitors.count_cells, selector)
        return size

    def select_icoords(self, SelectorObject selector, np.uint64_t size = -1):
        cdef GridVisitorData data
        self.setup_data(&data)
        if size == -1:
            size = 0
            data.array = <void*>(&size)
            self.visit_grids(&data,  grid_visitors.count_cells, selector)
        cdef np.ndarray[np.int64_t, ndim=2] icoords 
        icoords = np.zeros((size, 3), dtype="int64")
        data.array = icoords.data
        self.visit_grids(&data, grid_visitors.icoords_cells, selector)
        return icoords

    def select_ires(self, SelectorObject selector, np.uint64_t size = -1):
        cdef GridVisitorData data
        self.setup_data(&data)
        if size == -1:
            size = 0
            data.array = <void*>(&size)
            self.visit_grids(&data,  grid_visitors.count_cells, selector)
        cdef np.ndarray[np.int64_t, ndim=1] ires 
        ires = np.zeros(size, dtype="int64")
        data.array = ires.data
        self.visit_grids(&data, grid_visitors.ires_cells, selector)
        return ires

    def select_fcoords(self, SelectorObject selector, np.uint64_t size = -1):
        cdef GridVisitorData data
        self.setup_data(&data)
        if size == -1:
            size = 0
            data.array = <void*>(&size)
            self.visit_grids(&data,  grid_visitors.count_cells, selector)
        cdef np.ndarray[np.float64_t, ndim=2] fcoords 
        fcoords = np.zeros((size, 3), dtype="float64")
        data.array = fcoords.data
        self.visit_grids(&data, grid_visitors.fcoords_cells, selector)
        return fcoords

    def select_fwidth(self, SelectorObject selector, np.uint64_t size = -1):
        cdef GridVisitorData data
        self.setup_data(&data)
        if size == -1:
            size = 0
            data.array = <void*>(&size)
            self.visit_grids(&data,  grid_visitors.count_cells, selector)
        cdef np.ndarray[np.float64_t, ndim=2] fwidth 
        fwidth = np.zeros((size, 3), dtype="float64")
        data.array = fwidth.data
        self.visit_grids(&data, grid_visitors.fwidth_cells, selector)
        return fwidth
    
cdef class MatchPointsToGrids:

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __cinit__(self, GridTree tree,
                  int num_points, 
                  np.ndarray[np.float64_t, ndim=1] x,
                  np.ndarray[np.float64_t, ndim=1] y,
                  np.ndarray[np.float64_t, ndim=1] z):

        cdef int i
        
        self.num_points = num_points
        self.xp = <np.float64_t *> malloc(
                sizeof(np.float64_t) * num_points)
        self.yp = <np.float64_t *> malloc(
                sizeof(np.float64_t) * num_points)
        self.zp = <np.float64_t *> malloc(
                sizeof(np.float64_t) * num_points)
        self.point_grids = <np.int64_t *> malloc(
                sizeof(np.int64_t) * num_points)
        for i in range(num_points):
            self.xp[i] = x[i]
            self.yp[i] = y[i]
            self.zp[i] = z[i]
            self.point_grids[i] = -1
        self.tree = tree

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def find_points_in_tree(self):
        cdef np.ndarray[np.int64_t, ndim=1] pt_grids
        cdef int i, j
        cdef np.uint8_t in_grid
        pt_grids = np.zeros(self.num_points, dtype='int64')
        for i in range(self.num_points):
            in_grid = 0
            for j in range(self.tree.num_root_grids):
                if not in_grid: 
                    in_grid = self.check_position(i, self.xp[i], self.yp[i], self.zp[i],
                                                  &self.tree.root_grids[j])
        for i in range(self.num_points):
            pt_grids[i] = self.point_grids[i]
        return pt_grids

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef np.uint8_t check_position(self,
                                   np.int64_t pt_index, 
                                   np.float64_t x,
                                   np.float64_t y,
                                   np.float64_t z,
                                   GridTreeNode * grid):
        cdef int i
        cdef np.uint8_t in_grid
        in_grid = self.is_in_grid(x, y, z, grid)
        if in_grid:
            if grid.num_children > 0:
                in_grid = 0
                for i in range(grid.num_children):
                    if not in_grid:
                        in_grid = self.check_position(pt_index, x, y, z, grid.children[i])
                if not in_grid:
                    self.point_grids[pt_index] = grid.index
                    in_grid = 1
            else:
                self.point_grids[pt_index] = grid.index
                in_grid = 1
        return in_grid
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.uint8_t is_in_grid(self,
             np.float64_t x,
             np.float64_t y,
             np.float64_t z,
             GridTreeNode * grid):
        if x >= grid.right_edge[0]: return 0
        if y >= grid.right_edge[1]: return 0
        if z >= grid.right_edge[2]: return 0
        if x < grid.left_edge[0]: return 0
        if y < grid.left_edge[1]: return 0
        if z < grid.left_edge[2]: return 0
        return 1

