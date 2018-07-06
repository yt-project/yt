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
        The positions of particles in kdtree sorted order. Currently assumed
        to be 3D postions.
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
    cdef np.float64_t tpos, ma, sq_dist
    cdef np.float64_t* pos
    cdef np.float64_t[:] smoothing_length = np.empty(n_particles)
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
            process_node_points(leafnode, pos, input_positions, queue, i)

            # Traverse the rest of the kdtree to finish the neighbor list
            find_knn(
                c_tree.root, queue, input_positions, pos, leafnode.leafid, i)

            smoothing_length[i] = sqrt(queue.heap_ptr[0])

    pbar.update(n_particles-1)
    pbar.finish()
    return np.asarray(smoothing_length)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def generate_nn_list(np.int64_t[:, :, :, ::1] pids, np.float64_t[:, :, :, ::1] dists,
                     int offset, PyKDTree kdtree,
                     np.int64_t[:] tree_id, np.float64_t[:] bounds,
                     np.int64_t[:] size, np.float64_t[:,::1] input_pos):
    """Calculate array of distances to the nth nearest neighbor

    Parameters
    ----------

    Returns
    -------

    """
    cdef KDTree* c_tree = kdtree._tree
    cdef Node* leafnode
    cdef np.float64_t* pos
    cdef int i, j, k, p
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(pids.shape[3], True)
    cdef np.float64_t dx, dy, dz
    cdef np.float64_t[::1] voxel_position = np.zeros(shape=(3))

    dx = (bounds[1] - bounds[0]) / size[0]
    dy = (bounds[3] - bounds[2]) / size[1]
    dz = (bounds[5] - bounds[4]) / size[2]

    with nogil:
        for i in range(size[0]):
            for j in range(size[1]):
                for k in range(size[2]):
                    if i % CHUNKSIZE == 0:
                        with gil:
                            PyErr_CheckSignals()

                    voxel_position[0] = bounds[0] + (i+0.5)*dx
                    voxel_position[1] = bounds[2] + (j+0.5)*dy
                    voxel_position[2] = bounds[4] + (k+0.5)*dz

                    pos = &(voxel_position[0])
                    leafnode = c_tree.search(&pos[0])

                    queue.size = 0
                    queue.heap[:] = dists[i, j, k, :]
                    queue.pids[:] = pids[i, j, k, :]

                    if offset > 0:
                        queue.reset()
                        queue.size = queue.max_elements

                    process_node_points_pid(leafnode, pos, input_pos,
                                                offset, tree_id, queue)

                    # Traverse the rest of the kdtree to finish the neighbor
                    # list
                    find_knn_pid(c_tree.root, queue, input_pos, offset, tree_id,
                                 pos, leafnode.leafid)

                    dists[i, j, k, :] = queue.heap[:]
                    pids[i, j, k, :] = queue.pids[:]



@cython.boundscheck(False)
@cython.wraparound(False)
cdef int find_knn(Node* node,
                  BoundedPriorityQueue queue,
                  np.float64_t[:, ::1] positions,
                  np.float64_t* pos,
                  uint32_t skipleaf,
                  uint64_t skipidx) nogil except -1:
    # if we aren't a leaf then we keep travsersing until we find a leaf, else we
    # we actually begin to check the leaf
    if not node.is_leaf:
        if not cull_node(node.less, pos, queue, skipleaf):
            find_knn(node.less, queue, positions, pos, skipleaf, skipidx)
        if not cull_node(node.greater, pos, queue, skipleaf):
            find_knn(node.greater, queue, positions, pos, skipleaf, skipidx)
    else:
        if not cull_node(node, pos, queue, skipleaf):
            process_node_points(node, pos, positions, queue, skipidx)
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int find_knn_pid(Node* node,
                  BoundedPriorityQueue queue,
                  np.float64_t[:, ::1] positions,
                  np.int64_t offset,
                  np.int64_t[:] tree_id,
                  np.float64_t* pos,
                  uint32_t skipleaf) nogil except -1:

    # if we aren't a leaf then we keep travsersing until we find a leaf, else we
    # we actually begin to check the leaf
    if not node.is_leaf:
        if not cull_node_pid(node.less, pos, queue, skipleaf) and (tree_id[node.left_idx] - offset) > 0:
            find_knn_pid(node.less, queue, positions, offset, tree_id, pos, skipleaf)
        if not cull_node_pid(node.greater, pos, queue, skipleaf) and (tree_id[node.left_idx + node.children] - offset) < positions.shape[0]:
            find_knn_pid(node.greater, queue, positions, offset, tree_id, pos, skipleaf)
    else:
        if not cull_node_pid(node, pos, queue, skipleaf):
            process_node_points_pid(node, pos, positions, offset, tree_id, queue)
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int cull_node(Node* node,
                          np.float64_t* pos,
                          BoundedPriorityQueue queue,
                          uint32_t skipleaf) nogil except -1:
    # this function essentially checks if the node can possibly have particles
    # which are nearest neighbours, if it does, then it returns False. if it can
    # be skipped, it returns True
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
    return (ndist > queue.heap[0] and queue.size == queue.max_elements)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int cull_node_pid(Node* node,
                          np.float64_t* pos,
                          BoundedPriorityQueue queue,
                          uint32_t skipleaf) nogil except -1:
    # this function essentially checks if the node can possibly have particles
    # which are nearest neighbours, if it does, then it returns False. if it can
    # be skipped, it returns True
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
    return (ndist > queue.heap_ptr[0] and queue.size == queue.max_elements)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int process_node_points(Node* node,
                                    np.float64_t* pos,
                                    np.float64_t[:, ::1] positions,
                                    BoundedPriorityQueue queue,
                                    uint64_t skip_idx) nogil except -1:
    cdef uint64_t i, idx
    cdef np.float64_t tpos, sq_dist
    cdef int j
    cdef np.float64_t* p_ptr
    for i in range(node.left_idx, node.left_idx + node.children):
        p_ptr = &(positions[i, 0])
        if i != skip_idx:
            sq_dist = 0
            for j in range(3):
                tpos = p_ptr[j] - pos[j]
                sq_dist += tpos*tpos
            queue.add(sq_dist)
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int process_node_points_pid(Node* node,
                                        np.float64_t* pos,
                                        np.float64_t[:, ::1] positions,
                                        np.int64_t offset,
                                        np.int64_t[:] tree_id,
                                        BoundedPriorityQueue queue,
                                        ) nogil except -1:
    cdef uint64_t i,idx_offset
    cdef np.float64_t tpos, sq_dist
    cdef int j
    for i in range(node.left_idx, node.left_idx + node.children):
        idx_offset = tree_id[i] - offset
        if(idx_offset < 0 or idx_offset > positions.shape[0]):
            continue
        sq_dist = 0
        for j in range(3):
            tpos = positions[idx_offset, j] - pos[j]
            sq_dist += tpos*tpos
        queue.add_pid(sq_dist, i)
    return 0
