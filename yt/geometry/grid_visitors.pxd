"""
Grid visitor definitions file




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np

cdef struct GridTreeNode:
    np.int32_t num_children
    np.int32_t level
    np.int64_t index
    np.float64_t left_edge[3]
    np.float64_t right_edge[3]
    GridTreeNode **children
    np.int64_t start_index[3]
    np.int32_t dims[3]
    np.float64_t dds[3]

cdef struct GridTreeNodePadded:
    np.int32_t num_children
    np.int32_t level
    np.int64_t index
    np.float64_t left_edge_x
    np.float64_t left_edge_y
    np.float64_t left_edge_z
    np.float64_t right_edge_x
    np.float64_t right_edge_y
    np.float64_t right_edge_z
    np.int64_t children_pointers
    np.int64_t start_index_x
    np.int64_t start_index_y
    np.int64_t start_index_z
    np.int32_t dims_x
    np.int32_t dims_y
    np.int32_t dims_z
    np.float64_t dds_x
    np.float64_t dds_y
    np.float64_t dds_z

cdef struct GridVisitorData:
    GridTreeNode *grid
    np.uint64_t index
    np.uint64_t global_index
    np.int64_t pos[3]       # position in ints
    int n_tuples
    int **child_tuples # [N_child][6], where 0-1 are x_start, x_end, etc.
    void *array
    int ref_factor # This may change on a grid-by-grid basis
                   # It is the number of cells a child grid has per dimension
                   # in a cell of this grid.
                            
cdef void free_tuples(GridVisitorData *data) nogil
cdef void setup_tuples(GridVisitorData *data) nogil
cdef np.uint8_t check_child_masked(GridVisitorData *data) nogil

ctypedef void grid_visitor_function(GridVisitorData *data,
                                         np.uint8_t selected) nogil
# This is similar in spirit to the way oct visitor functions work.  However,
# there are a few important differences.  Because the grid objects are expected
# to be bigger, we don't need to pass them along -- we will not be recursively
# visiting.  So the GridVisitorData will be updated in between grids.
# Furthermore, we're only going to use them for a much smaller subset of
# operations.  All child mask evaluation is going to be conducted inside the
# outermost level of the visitor function, and visitor functions will receive
# information about whether they have been selected and whether they are
# covered by child cells.

cdef grid_visitor_function count_cells
cdef grid_visitor_function mask_cells
cdef grid_visitor_function icoords_cells
cdef grid_visitor_function ires_cells
cdef grid_visitor_function fcoords_cells
cdef grid_visitor_function fwidth_cells
