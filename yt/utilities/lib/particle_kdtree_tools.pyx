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
def generate_nn_list(np.float64_t[:] bounds, np.int64_t[:] dimensions,
                     PyKDTree kdtree, int n_neighbors):
    """Calculate an array of distances to the nearest n_neighbours and which
    particle is that distance away.

    Parameters
    ----------

    bounds: arrays of floats with shape (6)
        The bounds of the region to  pixelize / voxelize.
    dimensions: the number of pixels / voxels to divide into.
    kdtree: A PyKDTree instance
        A kdtree to do nearest neighbors searches with
    n_neighbors: The neighbor number to calculate the distance to

    Returns
    -------

    nearest_neighbours: tuple of arrays of floats and ints with shape ( float (n_pixels,
    n_neighbours), int (n_pixels, n_neighbours) )

    """
    cdef KDTree* c_tree = kdtree._tree
    cdef Node* leafnode
    cdef uint64_t idx
    cdef uint32_t skipid
    cdef np.float64_t tpos, ma, sq_dist
    cdef np.float64_t[:] pos
    cdef uint64_t neighbor_id
    cdef int i, j, k, l, skip
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(n_neighbors, True)
    cdef uint64_t xsize, ysize, zsize, n_pixels
    cdef np.float64_t x_min, x_max, y_min, y_max, z_min, z_max, dx, dy, dz

    # setting up the pixels to loop through
    n_pixels = xsize*ysize*zsize
    xsize, ysize, zsize = dimensions[0], dimensions[1], dimensions[2]

    x_min = bounds[0]
    x_max = bounds[1]
    y_min = bounds[2]
    y_max = bounds[3]
    z_min = bounds[4]
    z_max = bounds[5]

    dx = (x_max - x_min) / xsize
    dy = (y_max - y_min) / ysize
    dz = (z_max - z_min) / zsize

    pos = np.array([x_min+dx/2, y_min+dy/2, z_min+dz/2])

    pbar = get_pbar("Generate nearest neighbours", n_pixels)
    with nogil:
        for i in range(0, xsize):
            for j in range(0, ysize):
                for k in range(0, zsize):
                    # reset queue to "empty" state, doing it this way avoids
                    # needing to reallocate memory
                    queue.size = 0

                    if i % CHUNKSIZE == 0:
                        with gil:
                            pbar.update(i-1)
                            PyErr_CheckSignals()

                    leafnode = c_tree.search(&pos[0])
                    skipid = leafnode.leafid

                    # Fill queue with particles in the node containing the
                    # particle we're searching for
                    #process_node_points(leafnode, pos, input_positions, queue,
                    #                    i)

                    # Traverse the rest of the kdtree to finish the neighbor
                    # list
                    #find_knn(c_tree.root, queue, input_positions, pos,
                    #         leafnode.leafid, i)

                    pos[0] += dx
                    pos[1] += dy
                    pos[2] += dz

    pbar.update(n_pixels-1)
    pbar.finish()

    return

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int find_knn(Node* node,
                  BoundedPriorityQueue queue,
                  np.float64_t[:, :] positions,
                  np.float64_t* pos,
                  uint32_t skipleaf,
                  uint64_t skipidx) nogil except -1:
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

