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

cdef class Node

cdef class Node:
    cdef public Node left
    cdef public Node right
    cdef public Node parent
    cdef public int grid
    cdef public np.int64_t node_id
    cdef public np.int64_t node_ind
    cdef np.float64_t left_edge[3]
    cdef np.float64_t right_edge[3]
    cdef public data
    cdef Split * split
    cdef int level

cdef int point_in_node(Node node, np.ndarray[np.float64_t, ndim=1] point)
cdef Node _find_node(Node node, np.float64_t *point)
cdef int _kd_is_leaf(Node node)
