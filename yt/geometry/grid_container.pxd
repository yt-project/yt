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

from libc.stdlib cimport malloc, free
from libc.math cimport nearbyint, rint
from yt.geometry.selection_routines cimport SelectorObject, _ensure_code
from yt.utilities.lib.fp_utils cimport iclip

cdef struct GridTreeNode :
    int num_children
    int level
    int index
    np.float64_t left_edge[3]
    np.float64_t right_edge[3]
    GridTreeNode **children
                
cdef class GridTree:
    cdef GridTreeNode *grids
    cdef GridTreeNode *root_grids
    cdef int num_grids
    cdef int num_root_grids
    cdef int num_leaf_grids

cdef class MatchPointsToGrids :

    cdef int num_points
    cdef np.float64_t * xp
    cdef np.float64_t * yp
    cdef np.float64_t * zp
    cdef GridTree tree
    cdef np.int64_t * point_grids
    cdef np.uint8_t check_position(self,
                                   np.int64_t pt_index, 
                                   np.float64_t x,
                                   np.float64_t y,
                                   np.float64_t z,
                                   GridTreeNode * grid)

    cdef np.uint8_t is_in_grid(self,
			 np.float64_t x,
			 np.float64_t y,
			 np.float64_t z,
			 GridTreeNode * grid)

cdef class FastGridSelectionHelper:
    cdef public object index
    cdef public object grid_ind

