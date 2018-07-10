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
@cython.wraparound(True)
@cython.cdivision(True)
def generate_nn_list(np.int64_t[:, :, :, ::1] pids, np.float64_t[:, :, :, ::1] dists,
                     PyKDTree kdtree, np.float64_t[:] bounds, np.int64_t[:] size,
                     np.float64_t[:,::1] input_pos, int gather_type=0, int offset=0):
    """Calculate array of distances to the nth nearest neighbors, by recursively
    searching the tree

    Parameters
    ----------

    pids:   Array of ints to store the identity of the particle at a certain
            distance (voxels_x, voxels_y, voxels_z, n_neighbors)
    dists:  Array of floats to store the distance of a nearest neighbor distance
            (voxels_x, voxels_y, voxels_z, n_neighbors)
    kdtree: A PyKDTree instance
    bounds: The region for the voxels to cover in the format x_min, xmax, ymin,
            ... (6)
    size:   The number of voxels to have in each dimension
    input_pos: Array of floats with shape (n_particles, 3)
            The positions of particles in kdtree sorted order. Currently assumed
            to be 3D postions.
    gather_type: An integer with 0 meaning that all dimensions are used in the
            distance and 1 means only the first two are used. This would be 1
            in a slice plot.
    offset: This is the number of gas particle previously calculated and this is
            to identify a particle from the tree id

    Returns
    -------
    No returns, the function mutates the dists and the pids

    """
    cdef KDTree* c_tree = kdtree._tree
    cdef Node* leafnode
    cdef np.float64_t* pos
    cdef int i, j, k, p
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(pids.shape[3], True)
    cdef np.float64_t dx, dy, dz
    cdef np.float64_t[::1] voxel_position = np.zeros(shape=(3))
    cdef np.int64_t[:] tree_id

    tree_id = kdtree.idx.astype("int64")

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

                    queue.heap[:] = dists[i, j, k, :]
                    queue.pids[:] = pids[i, j, k, :]
                    queue.size = pids.shape[3]

                    pos = &(voxel_position[0])
                    leafnode = c_tree.search(&pos[0])

                    # Traverse the rest of the kdtree to finish the neighbor
                    # list
                    find_knn_pid(c_tree.root, queue, input_pos, offset, tree_id,
                                pos, leafnode.leafid, &(bounds[0]), gather_type)

                    dists[i, j, k, :] = queue.heap[:]
                    pids[i, j, k, :] = queue.pids[:]

@cython.boundscheck(False)
@cython.wraparound(True)
@cython.cdivision(True)
def generate_nn_list_guess(np.int64_t[:, :, :, ::1] pids, np.float64_t[:, :, :, ::1] dists,
                     np.int64_t[:, : , :] queue_sizes, PyKDTree kdtree,
                     np.float64_t[:] bounds, np.int64_t[:] size,
                     np.float64_t[:,::1] input_pos, int offset=0, int gather_type=0):
    """Calculate array of distances to the nth nearest neighbor from the
    leafnode the voxel centre is in, i.e, the best guess of neighbors

    Parameters
    ----------

    pids:   Array of ints to store the identity of the particle at a certain
            distance (voxels_x, voxels_y, voxels_z, n_neighbors)
    dists:  Array of floats to store the distance of a nearest neighbor distance
            (voxels_x, voxels_y, voxels_z, n_neighbors)
    queue_sizes: The current size of the queue used in the initial chunk loading
    kdtree: A PyKDTree instance
    bounds: The region for the voxels to cover in the format x_min, xmax, ymin,
            ... (6)
    size:   The number of voxels to have in each dimension
    input_pos: Array of floats with shape (n_particles, 3)
            The positions of particles in kdtree sorted order. Currently assumed
            to be 3D postions.
    gather_type: An integer with 0 meaning that all dimensions are used in the
            distance and 1 means only the first two are used
    offset: This is the number of gas particle previously calculated and this is
            to identify a particle from the tree id

    Returns
    -------
    No returns, the function mutates the dists, pids and queue_sizes

    """
    cdef KDTree* c_tree = kdtree._tree
    cdef Node* leafnode
    cdef np.float64_t* pos
    cdef int i, j, k, p
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(pids.shape[3], True)
    cdef np.float64_t dx, dy, dz
    cdef np.float64_t[::1] voxel_position = np.zeros(shape=(3))
    cdef np.int64_t[:] tree_id

    tree_id = kdtree.idx.astype("int64")

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

                    queue.heap[:] = dists[i, j, k, :]
                    queue.pids[:] = pids[i, j, k, :]
                    queue.size = queue_sizes[i, j, k]

                    pos = &(voxel_position[0])
                    leafnode = c_tree.search(&pos[0])
                    process_node_points_pid(leafnode, pos, input_pos,
                                offset, tree_id, queue, &(bounds[0]), gather_type)

                    queue_sizes[i, j, k] = queue.size
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
    # if we aren't a leaf then we keep traversing until we find a leaf, else we
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
                  uint32_t skipleaf,
                  np.float64_t * bounds,
                  int gather_type) nogil except -1:

    if not node.is_leaf:
        if not cull_node(node.less, pos, queue, skipleaf, gather_type):
            find_knn_pid(node.less, queue, positions, offset, tree_id, pos,
                         skipleaf, bounds, gather_type)
        if not cull_node(node.greater, pos, queue, skipleaf, gather_type):
            find_knn_pid(node.greater, queue, positions, offset, tree_id, pos,
                         skipleaf, bounds, gather_type)
    else:
        if not cull_node(node, pos, queue, skipleaf, gather_type):
            process_node_points_pid(node, pos, positions, offset, tree_id,
                                    queue, bounds, gather_type)
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int cull_node(Node* node,
                          np.float64_t* pos,
                          BoundedPriorityQueue queue,
                          uint32_t skipleaf,
                          int gather_type=0) nogil except -1:
    cdef int dim = 3
    if gather_type == 1:
        dim = 2

    cdef np.float64_t v
    cdef np.float64_t tpos, ndist = 0
    cdef uint32_t leafid
    if node.leafid == skipleaf:
        return True
    for k in range(dim):
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
                                        np.float64_t * bounds,
                                        int gather_type) nogil except -1:
    # if gather type is 1, then we ignore the third coordinate for a slice
    cdef int dim = 3
    if gather_type == 1:
        dim = 2

    cdef uint64_t i,idx_offset
    cdef np.float64_t tpos, sq_dist
    cdef int j
    for i in range(node.left_idx, node.left_idx + node.children):
        idx_offset = tree_id[i] - offset
        if(idx_offset < 0 or idx_offset > positions.shape[0]):
            continue

        # skip the particle if it is not in the region
        if positions[idx_offset, 0] < bounds[0] or \
            positions[idx_offset, 0] > bounds[1]:
                continue
        if positions[idx_offset, 1] < bounds[2] or \
            positions[idx_offset, 1] > bounds[3]:
                continue
        if positions[idx_offset, 2] < bounds[4] or \
            positions[idx_offset, 2] > bounds[5]:
                continue

        sq_dist = 0
        for j in range(dim):
            tpos = positions[idx_offset, j] - pos[j]
            sq_dist += tpos*tpos

        queue.add_pid(sq_dist, i)

    return 0
