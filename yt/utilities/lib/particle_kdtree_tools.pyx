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
                              PyKDTree kdtree, int n_neighbors
                              ):
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

    # these are things which are not necessary for this function, but we set
    # them up to pass to the utility functions as we use them elsewhere
    cdef int offset = 0
    cdef np.int64_t[:] tree_id = np.zeros((1), dtype="int64") - 1

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
            process_node_points(leafnode, queue, input_positions, pos, i,
                                tree_id, offset)

            # Traverse the rest of the kdtree to finish the neighbor list
            find_knn(c_tree.root, queue, input_positions, pos, leafnode.leafid,
                     i, tree_id, offset)

            smoothing_length[i] = sqrt(queue.heap_ptr[0])

    pbar.update(n_particles-1)
    pbar.finish()
    return np.asarray(smoothing_length)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def query(PyKDTree kdtree, np.float64_t [:, ::1] tree_positions,
          np.float64_t[:, ::1] input_positions, int num_neigh
          ):
    """This is a KD nearest neighbor search with a similar interface to the
    scikdtree used for performance comparisons

    Parameters
    ----------

    kdtree: A PyKDTree instance
        A kdtree to do nearest neighbors searches with
    tree_positions: arrays of floats with shape (n_particles, 3)
        The positions of particles in kdtree non-sorted order. Currently assumed
        to be 3D postions.
    input_positions: array of floats (any, 3)
        The positions to gather the nearest neighbors.
    kdtree: A PyKDTree instance
        A kdtree to do nearest neighbors searches with
    num_neigh: The neighbor number to calculate the distance to

    Returns
    -------

    dists: arrays of floats with shape (n_particles, num_neigh)
        The the nearest neighbor distances
    pids: arrays of ints with shape (n_particles, num_neigh)
        The particle ids of the nearest neighbors
    """

    cdef KDTree* c_tree = kdtree._tree
    cdef Node* leafnode
    cdef np.float64_t* pos
    cdef int i, skipidx = -1, offset = 0
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(num_neigh, True)
    cdef np.int64_t[:] tree_id
    cdef np.float64_t[:, :] dists
    cdef np.int64_t[:, :] pids

    dists = np.zeros((input_positions.shape[0], num_neigh), dtype="float64")
    pids = np.zeros((input_positions.shape[0], num_neigh), dtype="int64")

    # if we haven't sorted then we use the tree ids
    tree_id = kdtree.idx.astype("int64")

    for i in range(0, input_positions.shape[0]):
        queue.size = 0

        pos = &(input_positions[i, 0])
        leafnode = c_tree.search(pos)

        process_node_points(leafnode, queue, tree_positions, pos, skipidx, tree_id,
                            offset)
        find_knn(c_tree.root, queue, tree_positions, pos, leafnode.leafid,
                 skipidx, tree_id, offset)

        dists[i, :] = queue.heap[:]
        pids[i, :] = queue.pids[:]

    return np.asarray(dists), np.asarray(pids)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def knn_list_guess(np.float64_t [:, ::1] tree_positions, np.float64_t[:, :, :, ::1] dists,
                   np.int64_t[:, :, :, ::1] pids, np.int64_t[:, :, :] q_sizes,
                   PyKDTree kdtree, np.float64_t[:] bounds, np.int64_t[:] size,
                   int num_neigh, int sorted=0, int offset=0
                   ):
    """This is a KD nearest neighbor search to calculate the first guess, this
    is useful for filling a nearest neighbor list chunkwise for a voxel.

    Parameters
    ----------

    tree_positions: arrays of floats with shape (n_particles, 3)
        The positions of particles in kdtree non-sorted order. Currently assumed
        to be 3D postions.
    kdtree: A PyKDTree instance
        A kdtree to do nearest neighbors searches with
    dists: arrays of floats with shape (n_particles, num_neigh)
        The the nearest neighbor distances
    pids: arrays of ints with shape (n_particles, num_neigh)
        The particle ids of the nearest neighbors
    kdtree: A PyKDTree instance
        A kdtree to do nearest neighbors searches with
    bounds: array of floats
        The boundaries of the grid.
    size: array of ints
        The number of cells to divide the grid into
    num_neigh: The neighbor number to calculate the distance to
    sorted: int
        have the tree input positions been sorted with the tree idx?
    offset: int
        the number of particles already been filled, this is used for chunkwise
        neighbor finding
    """

    # the positions are the positions to find the k nn and the dists and pids
    cdef KDTree* c_tree = kdtree._tree
    cdef Node* leafnode
    cdef np.float64_t* pos
    cdef int i, j, k, skipidx = -1
    cdef double dx, dy, dz
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(num_neigh, True)
    cdef np.int64_t[:] tree_id
    cdef np.float64_t[:] voxel_pos = np.zeros(3, dtype="float64")

    dx = (bounds[1] - bounds[0]) / size[0]
    dy = (bounds[3] - bounds[2]) / size[1]
    dz = (bounds[5] - bounds[4]) / size[2]

    # if we haven't sorted then we use the tree ids to map the particle ids
    if not sorted:
        tree_id = kdtree.idx.astype("int64")
    else:
        tree_id = np.zeros((1), dtype="int64") - 1

    for i in range(0, size[0]):
        for j in range(0, size[1]):
            for k in range(0, size[2]):
                queue.size = q_sizes[i, j, k]
                queue.heap[:] = dists[i, j, k, :]
                queue.pids[:] = pids[i, j, k, :]

                voxel_pos[0] = bounds[0] + (i+0.5)*dx
                voxel_pos[1] = bounds[2] + (j+0.5)*dy
                voxel_pos[2] = bounds[4] + (k+0.5)*dz

                pos = &(voxel_pos[0])
                leafnode = c_tree.search(pos)

                process_node_points(leafnode, queue, tree_positions, pos, skipidx,
                            tree_id, offset)

                q_sizes[i, j, k] = queue.size
                dists[i, j, k, :] = queue.heap[:]
                pids[i, j, k, :] = queue.pids[:]

    return

@cython.boundscheck(False)
@cython.wraparound(True)
@cython.cdivision(True)
def knn_list(np.float64_t [:, ::1] tree_positions, np.float64_t[:, :, :, ::1] dists,
             np.int64_t[:, :, :, ::1] pids,  PyKDTree kdtree, np.float64_t[:] bounds,
             np.int64_t[:] size, int num_neigh, int sorted=0, int offset=0
             ):
    """This is a KD nearest neighbor search to calculate improve the first guess
    with tree traversing, this is useful for filling a nearest neighbor list
    chunkwise

    Parameters
    ----------

    tree_positions: arrays of floats with shape (n_particles, 3)
        The positions of particles in kdtree non-sorted order. Currently assumed
        to be 3D postions.
    kdtree: A PyKDTree instance
        A kdtree to do nearest neighbors searches with
    dists: arrays of floats with shape (n_particles, num_neigh)
        The the nearest neighbor distances
    pids: arrays of ints with shape (n_particles, num_neigh)
        The particle ids of the nearest neighbors
    kdtree: A PyKDTree instance
        A kdtree to do nearest neighbors searches with
    input_positions: array of floats (any, 3)
        The positions to gather the nearest neighbors.
    num_neigh: The neighbor number to calculate the distance to
    sorted: int
        have the tree input positions been sorted with the tree idx?
    offset: int
        the number of particles already been filled, this is used for chunkwise
        neighbor finding
    """

    # the positions are the positions to find the k nn and the dists and pids
    cdef KDTree* c_tree = kdtree._tree
    cdef Node* leafnode
    cdef np.float64_t* pos
    cdef int i, j, k, skipidx = -1
    cdef double dx, dy, dz
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(num_neigh, True)
    cdef np.int64_t[:] tree_id
    cdef np.float64_t[:] voxel_pos = np.zeros(3, dtype="float64")

    dx = (bounds[1] - bounds[0]) / size[0]
    dy = (bounds[3] - bounds[2]) / size[1]
    dz = (bounds[5] - bounds[4]) / size[2]

    # if we haven't sorted then we use the tree ids to map the particle ids
    if not sorted:
        tree_id = kdtree.idx.astype("int64")
    else:
        tree_id = np.zeros((1), dtype="int64") - 1

    for i in range(0, size[0]):
        for j in range(0, size[1]):
            for k in range(0, size[2]):
                queue.size = queue.max_elements
                queue.heap[:] = dists[i, j, k, :]
                queue.pids[:] = pids[i, j, k, :]

                voxel_pos[0] = bounds[0] + (i+0.5)*dx
                voxel_pos[1] = bounds[2] + (j+0.5)*dy
                voxel_pos[2] = bounds[4] + (k+0.5)*dz

                pos = &(voxel_pos[0])
                leafnode = c_tree.search(pos)

                find_knn(c_tree.root, queue, tree_positions, pos, leafnode.leafid,
                         skipidx, tree_id, offset)

                dists[i, j, k, :] = queue.heap[:]
                pids[i, j, k, :] = queue.pids[:]

    return

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int find_knn(Node* node,
                  BoundedPriorityQueue queue,
                  np.float64_t[:, ::1] tree_positions,
                  np.float64_t* pos,
                  uint32_t skipleaf,
                  uint64_t skipidx,
                  np.int64_t[:] tree_id,
                  int offset
                  ) nogil except -1:
    # if we aren't a leaf then we keep traversing until we find a leaf, else we
    # we actually begin to check the leaf
    if not node.is_leaf:
        if not cull_node(node.less, pos, queue, skipleaf):
            find_knn(node.less, queue, tree_positions, pos, skipleaf, skipidx,
                     tree_id, offset)
        if not cull_node(node.greater, pos, queue, skipleaf):
            find_knn(node.greater, queue, tree_positions, pos, skipleaf, skipidx,
                     tree_id, offset)
    else:
        if not cull_node(node, pos, queue, skipleaf):
            process_node_points(node, queue, tree_positions, pos, skipidx, tree_id,
                                offset)
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int cull_node(Node* node,
                          np.float64_t* pos,
                          BoundedPriorityQueue queue,
                          uint32_t skipleaf,
                          ) nogil except -1:
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
cdef inline int process_node_points(Node* node,
                                    BoundedPriorityQueue queue,
                                    np.float64_t[:, ::1] positions,
                                    np.float64_t* pos,
                                    int skipidx,
                                    np.int64_t[:] tree_id,
                                    int offset,
                                    ) nogil except -1:
    cdef uint64_t i, idx_offset
    cdef np.float64_t tpos, sq_dist
    cdef int j

    for i in range(node.left_idx, node.left_idx + node.children):
        if tree_id[0] == -1:
            idx_offset = i
        else:
            idx_offset = tree_id[i] - offset

        if(idx_offset < 0 or idx_offset > positions.shape[0] or \
           idx_offset == skipidx):
            continue

        sq_dist = 0
        for j in range(3):
            tpos = positions[idx_offset, j] - pos[j]
            sq_dist += tpos*tpos

        queue.add_pid(sq_dist, i)
    return 0

