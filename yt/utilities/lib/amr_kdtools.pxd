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

cimport numpy as np

cdef struct Split:
    int dim
    np.float64_t pos

cdef class Node:
    cdef public Node left
    cdef public Node right
    cdef public Node parent
    cdef public int grid
    cdef public bint dirty
    cdef public np.int64_t node_id
    cdef public np.int64_t node_ind
    cdef np.float64_t left_edge[3]
    cdef np.float64_t right_edge[3]
    cdef public data
    cdef Split * split
    cdef int level
    cdef int point_in_node(self, np.float64_t[:] point)
    cdef Node _find_node(self, np.float64_t[:] point)
    cdef int _kd_is_leaf(self)
    cdef add_grid(self, np.float64_t[:] gle, np.float64_t[:] gre,
                       int gid,
                       int rank,
                       int size)
    cdef insert_grid(self,
                       np.float64_t[:] gle,
                       np.float64_t[:] gre,
                       int grid_id,
                       int rank,
                       int size)
    cpdef add_grids(self,
                       int ngrids,
                       np.float64_t[:,:] gles,
                       np.float64_t[:,:] gres,
                       np.int64_t[:] gids,
                       int rank,
                       int size)
    cdef int should_i_split(self, int rank, int size)
    cdef void insert_grids(self,
                       int ngrids,
                       np.float64_t[:,:] gles,
                       np.float64_t[:,:] gres,
                       np.int64_t[:] gids,
                       int rank,
                       int size)
    cdef split_grid(self,
                       np.float64_t[:] gle,
                       np.float64_t[:] gre,
                       int gid,
                       int rank,
                       int size)
    cdef int split_grids(self,
                       int ngrids,
                       np.float64_t[:,:] gles,
                       np.float64_t[:,:] gres,
                       np.int64_t[:] gids,
                       int rank,
                       int size)
    cdef geo_split(self,
                       np.float64_t[:] gle,
                       np.float64_t[:] gre,
                       int grid_id,
                       int rank,
                       int size)
    cdef void divide(self, Split * split)
