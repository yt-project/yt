"""
Cython tools for working with the PyKDTree particle KDTree.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np

cimport cython

from cpython.exc cimport PyErr_CheckSignals
from cykdtree.kdtree cimport PyKDTree, KDTree, Node, uint64_t, uint32_t

from libc.math cimport sqrt
from libcpp.vector cimport vector

from yt.funcs import get_pbar
from yt.utilities.lib.bounded_priority_queue cimport BoundedPriorityQueue

cdef int CHUNKSIZE = 4096

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def generate_smoothing_length(np.float64_t[:, ::1] input_positions,
                              PyKDTree kdtree,
                              int n_neighbors):
    """Calculate array of distances to the nth nearest neighbor

    Parameters
    ----------

    input_positions: arrays of floats with shape (n_particles, 3)
        The positions of particles. Current assumed to be 3D postions.
    kdtree: A PyKDTree instance
        A kdtree to do nearest neighbors searches with
    n_neighbors: The neighbor number to calculate the distance to

    Returns
    -------

    smoothing_lengths: arrays of flots with shape (n_particles, )
        The calculated smoothing lengths

    """
    cdef KDTree* c_tree = kdtree._tree
    cdef Node* leafnode
    cdef uint64_t idx
    cdef uint32_t skipid
    cdef int n_particles = input_positions.shape[0]
    cdef np.float64_t[:] smoothing_length = np.empty(n_particles)
    cdef np.float64_t tpos, ma, sq_dist
    cdef np.float64_t* pos
    cdef uint64_t neighbor_id
    cdef int i, j, k, l, skip
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(n_neighbors)

    pbar = get_pbar("Generate smoothing length", n_particles)
    with nogil:
        for i in range(n_particles):
            # reset queue to "empty" state, doing it this way avoids
            # needing to reallocate memory
            queue.size = 0

            if i % CHUNKSIZE == 0:
                with gil:
                    pbar.update(i-1)
                    PyErr_CheckSignals()

            pos = &(input_positions[i, 0])
            leafnode = c_tree.search(&pos[0])
            skipid = leafnode.leafid

            # Fill queue with particles in the node containing the particle
            # we're searching for
            process_node_points(leafnode, pos, input_positions, c_tree.all_idx,
                                queue, i)

            # Traverse the rest of the kdtree to finish the neighbor list
            find_knn(c_tree.root, queue, input_positions, pos, c_tree.all_idx,
                     leafnode.leafid, i)

            smoothing_length[i] = sqrt(queue.heap_ptr[0])

    pbar.update(n_particles-1)
    pbar.finish()
    return np.asarray(smoothing_length)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int find_knn(Node* node,
                  BoundedPriorityQueue queue,
                  np.float64_t[:, :] positions,
                  np.float64_t* pos,
                  uint64_t* all_idx,
                  uint32_t skipleaf,
                  uint64_t skipidx) nogil except -1:
    if not node.is_leaf:
        if not cull_node(node.less, pos, queue, skipleaf):
            find_knn(node.less, queue, positions, pos, all_idx, skipleaf,
                     skipidx)
        if not cull_node(node.greater, pos, queue, skipleaf):
            find_knn(node.greater, queue, positions, pos, all_idx, skipleaf,
                     skipidx)
    else:
        if not cull_node(node, pos, queue, skipleaf):
            process_node_points(node, pos, positions, all_idx, queue, skipidx)
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int cull_node(Node* node,
                          np.float64_t* pos,
                          BoundedPriorityQueue queue,
                          uint32_t skipleaf) nogil except -1:
    cdef np.float64_t v
    cdef np.float64_t tpos, ndist = 0
    cdef uint32_t leafid
    if node.leafid == skipleaf:
        return True
    for k in range(3):
        v = pos[k]
        if v < node.left_edge[k]:
            tpos = node.left_edge[k] - v
        elif v > node.right_edge[k]:
            tpos = v - node.right_edge[k]
        else:
            tpos = 0
        ndist += tpos*tpos
    return ndist > queue.heap_ptr[0]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int process_node_points(Node* node,
                                    np.float64_t* pos,
                                    np.float64_t[:, :] positions,
                                    uint64_t* all_idx,
                                    BoundedPriorityQueue queue,
                                    uint64_t skip_idx) nogil except -1:
    cdef uint64_t i, idx
    cdef np.float64_t tpos, sq_dist
    cdef int j
    cdef np.float64_t* p_ptr
    for i in range(node.left_idx, node.left_idx + node.children):
        idx = all_idx[i]
        p_ptr = &(positions[idx, 0])
        if idx != skip_idx:
            sq_dist = 0
            for j in range(3):
                tpos = p_ptr[j] - pos[j]
                sq_dist += tpos*tpos
            queue.add(sq_dist)
    return 0

