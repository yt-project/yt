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

# This structure allows the nearest neighbor finding to consider a subset of
# spatial dimensions, i.e the the spatial separation in the x and z coordinates
# could be consider by using set_axes_range(axes, 1), this would cause the while
# loops to skip the y dimensions, without the performance hit of an if statement
cdef struct axes_range:
    int start
    int stop
    int step

# skipaxis: x=0, y=1, z=2
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int set_axes_range(axes_range *axes, int skipaxis):
    axes.start = 0
    axes.stop = 3
    axes.step = 1
    if skipaxis == 0:
        axes.start = 1
    if skipaxis == 1:
        axes.step = 2
    if skipaxis == 2:
        axes.stop = 2
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def generate_smoothing_length(np.float64_t[:, ::1] input_positions,
                              PyKDTree kdtree, int n_neighbors):
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
    cdef np.int64_t skipaxis = -1

    cdef axes_range axes
    set_axes_range(&axes, -1)

    # these are things which are not necessary for this function, but we set
    # them up to pass to the utility functions as we use them elsewhere
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
                                &axes)

            # Traverse the rest of the kdtree to finish the neighbor list
            find_knn(c_tree.root, queue, input_positions, pos, leafnode.leafid,
                     i, &axes)

            smoothing_length[i] = sqrt(queue.heap_ptr[0])

    pbar.update(n_particles-1)
    pbar.finish()
    return np.asarray(smoothing_length)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def knn_position(np.float64_t [::1] position,
                      np.float64_t [:, ::1] tree_positions,
                      BoundedPriorityQueue queue,
                      PyKDTree kdtree,
                      np.int64_t skipaxis=-1):

    """This calcualtes the K nearest neighbors of an individual position

    Parameters
    ----------

    position: array of floats (3)
        The position to find the nearest neighbors
    tree_positions: arrays of floats with shape (n_particles, 3)
        The positions of particles in kdtree non-sorted order. Currently assumed
        to be 3D postions.
    queue: A BoundedPriorityQueue instance
        This prevents the costant reallocation
    kdtree: A PyKDTree instance
        A kdtree to do nearest neighbors searches with
    skipaxis: int
        Any physics dimensions which should be ignored when calculating
        distances, i.e in a projection plot
    """

    # the positions are the positions to find the k nn and the dists and pids
    cdef KDTree* c_tree = kdtree._tree
    cdef Node* leafnode
    cdef np.float64_t* pos
    cdef int skipidx = -1
    cdef int extra_nodes

    cdef axes_range axes
    set_axes_range(&axes, skipaxis)

    with nogil:
        queue.size = 0
        pos = &(position[0])

        leafnode = c_tree.search(pos)
        process_node_points(leafnode, queue, tree_positions, pos, skipidx, &axes)

        # if we want more neighbors than particles in a leaf, then we
        # loop through a few neighbors
        extra_nodes = 0
        while queue.size < queue.max_elements and extra_nodes < \
                leafnode.all_neighbors.size():
            process_node_points(c_tree.leaves[leafnode.all_neighbors[extra_nodes]],
                            queue, tree_positions, pos, skipidx, &axes)
            extra_nodes += 1

        find_knn(c_tree.root, queue, tree_positions, pos, leafnode.leafid, skipidx,
             &axes)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def knn_list_grid(np.float64_t [:, ::1] tree_positions, np.float64_t[:, :, :, ::1] dists,
             np.int64_t[:, :, :, ::1] pids,  PyKDTree kdtree, np.float64_t[:] bounds,
             np.int64_t[:] size, int num_neigh,
             np.int64_t skipaxis=-1):

    """This calcualtes the K nearest neighbors of a uniform grid with a bounds
    and size that have been inputted. This is useful for slice plots, arbitrary
    grids and projections.

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
    bounds: array of floats (6)
        The boundaries of the grid
    size: array of ints (3)
        The number of voxels to divide the grid into for each dimenions
    num_neigh: The neighbor number to calculate the distance to
    skipaxis: int
        Any physics dimensions which should be ignored when calculating
        distances, i.e in a projection plot
    """

    # the positions are the positions to find the k nn and the dists and pids
    cdef KDTree* c_tree = kdtree._tree
    cdef Node* leafnode
    cdef np.float64_t* pos
    cdef int i, j, k, p, q, skipidx = -1
    cdef double dx, dy, dz
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(num_neigh, True)
    cdef np.float64_t[:] voxel_pos = np.zeros(3, dtype="float64")
    cdef int extra_nodes

    cdef axes_range axes
    set_axes_range(&axes, skipaxis)

    dx = (bounds[1] - bounds[0]) / size[0]
    dy = (bounds[3] - bounds[2]) / size[1]
    dz = (bounds[5] - bounds[4]) / size[2]

    pbar = get_pbar("Generating neighbor lists", size[2]*size[0]*size[1])
    p = 0
    for i in range(0, size[0]):
        for j in range(0, size[1]):
            for k in range(0, size[2]):
                queue.size = 0

                voxel_pos[0] = bounds[0] + (i+0.5)*dx
                voxel_pos[1] = bounds[2] + (j+0.5)*dy
                voxel_pos[2] = bounds[4] + (k+0.5)*dz
                pos = &(voxel_pos[0])

                leafnode = c_tree.search(pos)
                process_node_points(leafnode, queue, tree_positions, pos,
                                    skipidx, &axes)

                # if we want more neighbors than particles in a leaf, then we
                # loop through a few neighbors
                extra_nodes = 0
                while queue.size < queue.max_elements and extra_nodes < \
                  leafnode.all_neighbors.size():
                    process_node_points(c_tree.leaves[leafnode.all_neighbors[extra_nodes]],
                                        queue, tree_positions, pos, skipidx, &axes)
                    extra_nodes += 1

                find_knn(c_tree.root, queue, tree_positions, pos, leafnode.leafid,
                         skipidx, &axes)

                dists[i, j, k, :] = queue.heap[:]
                pids[i, j, k, :] = queue.pids[:]

                p += 1
                if p % CHUNKSIZE:
                    pbar.update(CHUNKSIZE)
    pbar.finish()

    return

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int find_knn(Node* node,
                  BoundedPriorityQueue queue,
                  np.float64_t[:, ::1] tree_positions,
                  np.float64_t* pos,
                  uint32_t skipleaf,
                  uint64_t skipidx,
                  axes_range * axes,
                  ) nogil except -1:
    # if we aren't a leaf then we keep traversing until we find a leaf, else we
    # we actually begin to check the leaf
    if not node.is_leaf:
        if not cull_node(node.less, pos, queue, skipleaf, axes):
            find_knn(node.less, queue, tree_positions, pos, skipleaf, skipidx,
                     axes)
        if not cull_node(node.greater, pos, queue, skipleaf, axes):
            find_knn(node.greater, queue, tree_positions, pos, skipleaf,
                     skipidx, axes)
    else:
        if not cull_node(node, pos, queue, skipleaf, axes):
            process_node_points(node, queue, tree_positions, pos, skipidx,
                                axes)
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int cull_node(Node* node,
                          np.float64_t* pos,
                          BoundedPriorityQueue queue,
                          uint32_t skipleaf,
                          axes_range * axes,
                          ) nogil except -1:
    cdef int k
    cdef np.float64_t v
    cdef np.float64_t tpos, ndist = 0
    cdef uint32_t leafid

    if node.leafid == skipleaf:
        return True

    k = axes.start
    while k < axes.stop:
        v = pos[k]
        if v < node.left_edge[k]:
            tpos = node.left_edge[k] - v
        elif v > node.right_edge[k]:
            tpos = v - node.right_edge[k]
        else:
            tpos = 0
        ndist += tpos*tpos
        k += axes.step

    return (ndist > queue.heap[0] and queue.size == queue.max_elements)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int process_node_points(Node* node,
                                    BoundedPriorityQueue queue,
                                    np.float64_t[:, ::1] positions,
                                    np.float64_t* pos,
                                    int skipidx,
                                    axes_range * axes,
                                    ) nogil except -1:
    cdef uint64_t i, k
    cdef np.float64_t tpos, sq_dist
    for i in range(node.left_idx, node.left_idx + node.children):
        sq_dist = 0
        k = axes.start
        while k < axes.stop:
            tpos = positions[i, k] - pos[k]
            sq_dist += tpos*tpos
            k += axes.step

        queue.add_pid(sq_dist, i)
    return 0

