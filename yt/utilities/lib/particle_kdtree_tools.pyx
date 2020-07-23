# distutils: language = c++
"""
Cython tools for working with the PyKDTree particle KDTree.



"""


import numpy as np

cimport cython
cimport numpy as np
from cpython.exc cimport PyErr_CheckSignals
from libc.math cimport sqrt
from libcpp.vector cimport vector

from yt.utilities.lib.cykdtree.kdtree cimport KDTree, Node, PyKDTree, uint32_t, uint64_t

from yt.funcs import get_pbar

from yt.geometry.particle_deposit cimport get_kernel_func, kernel_func
from yt.utilities.lib.bounded_priority_queue cimport BoundedPriorityQueue, NeighborList


cdef int CHUNKSIZE = 4096

# This structure allows the nearest neighbor finding to consider a subset of
# spatial dimensions, i.e the spatial separation in the x and z coordinates
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
def generate_smoothing_length(np.float64_t[:, ::1] tree_positions,
                              PyKDTree kdtree, int n_neighbors):
    """Calculate array of distances to the nth nearest neighbor

    Parameters
    ----------

    tree_positions: arrays of floats with shape (n_particles, 3)
        The positions of particles in kdtree sorted order. Currently assumed
        to be 3D postions.
    kdtree: A PyKDTree instance
        A kdtree to do nearest neighbors searches with
    n_neighbors: The neighbor number to calculate the distance to

    Returns
    -------

    smoothing_lengths: arrays of floats with shape (n_particles, )
        The calculated smoothing lengths

    """
    cdef int i
    cdef KDTree * c_tree = kdtree._tree
    cdef int n_particles = tree_positions.shape[0]
    cdef np.float64_t * pos
    cdef np.float64_t[:] smoothing_length = np.empty(n_particles)
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(n_neighbors)
    cdef np.int64_t skipaxis = -1

    # We are using all spatial dimensions
    cdef axes_range axes
    set_axes_range(&axes, -1)

    pbar = get_pbar("Generate smoothing length", n_particles)
    with nogil:
        for i in range(n_particles):
            # Reset queue to "empty" state, doing it this way avoids
            # needing to reallocate memory
            queue.size = 0

            if i % CHUNKSIZE == 0:
                with gil:
                    pbar.update(i-1)
                    PyErr_CheckSignals()

            pos = &(tree_positions[i, 0])
            find_neighbors(pos, tree_positions, queue, c_tree, i, &axes)

            smoothing_length[i] = sqrt(queue.heap_ptr[0])

    pbar.update(n_particles-1)
    pbar.finish()
    return np.asarray(smoothing_length)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def estimate_density(np.float64_t[:, ::1] tree_positions, np.float64_t[:] mass,
                      np.float64_t[:] smoothing_length,
                      PyKDTree kdtree, kernel_name="cubic"):
    """Estimate density using SPH gather method.

    Parameters
    ----------

    tree_positions: array of floats with shape (n_particles, 3)
        The positions of particles in kdtree sorted order. Currently assumed
        to be 3D postions.
    mass: array of floats with shape (n_particles)
        The masses of particles in kdtree sorted order.
    smoothing_length: array of floats with shape (n_particles)
        The smoothing lengths of particles in kdtree sorted order.
    kdtree: A PyKDTree instance
        A kdtree to do nearest neighbors searches with.
    kernel_name: str
        The name of the kernel function to use in density estimation.

    Returns
    -------

    density: array of floats with shape (n_particles)
        The calculated density.

    """
    cdef int i, j, k
    cdef KDTree * c_tree = kdtree._tree
    cdef int n_particles = tree_positions.shape[0]
    cdef np.float64_t h_i2, ih_i2, q_ij
    cdef np.float64_t * pos
    cdef np.float64_t[:] density = np.empty(n_particles)
    cdef kernel_func kernel = get_kernel_func(kernel_name)
    cdef NeighborList nblist = NeighborList()

    # We are using all spatial dimensions
    cdef axes_range axes
    set_axes_range(&axes, -1)

    pbar = get_pbar("Estimating density", n_particles)
    with nogil:
        for i in range(n_particles):
            # Reset list to "empty" state, doing it this way avoids
            # needing to reallocate memory
            nblist.size = 0

            if i % CHUNKSIZE == 0:
                with gil:
                    pbar.update(i - 1)
                    PyErr_CheckSignals()

            pos = &(tree_positions[i, 0])
            h_i2 = smoothing_length[i] ** 2
            find_neighbors_ball(pos, h_i2, tree_positions, nblist, c_tree, i, &axes)
            ih_i2 = 1.0 / h_i2

            # See eq. 10 of Price 2012
            density[i] = mass[i] * kernel(0)
            for k in range(nblist.size):
                j = nblist.pids[k]
                q_ij = sqrt(nblist.data[k] * ih_i2)
                density[i] += mass[j] * kernel(q_ij)

    pbar.update(n_particles - 1)
    pbar.finish()
    return np.asarray(density)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int find_neighbors(np.float64_t * pos, np.float64_t[:, ::1] tree_positions,
                        BoundedPriorityQueue queue, KDTree * c_tree,
                        uint64_t skipidx, axes_range * axes) nogil except -1:
    cdef Node* leafnode

    # Make an initial guess based on the closest node
    leafnode = c_tree.search(&pos[0])
    process_node_points(leafnode, queue, tree_positions, pos, skipidx, axes)

    # Traverse the rest of the kdtree to finish the neighbor list
    find_knn(c_tree.root, queue, tree_positions, pos, leafnode.leafid, skipidx,
             axes)

    return 0

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
        if i == skipidx:
            continue

        sq_dist = 0.0

        k = axes.start
        while k < axes.stop:
            tpos = positions[i, k] - pos[k]
            sq_dist += tpos*tpos
            k += axes.step

        queue.add_pid(sq_dist, i)

    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int find_neighbors_ball(np.float64_t * pos, np.float64_t r2,
                             np.float64_t[:, ::1] tree_positions,
                             NeighborList nblist, KDTree * c_tree,
                             uint64_t skipidx, axes_range * axes
                             ) nogil except -1:
    """Find neighbors within a ball."""
    cdef Node* leafnode

    # Make an initial guess based on the closest node
    leafnode = c_tree.search(&pos[0])
    process_node_points_ball(leafnode, nblist, tree_positions, pos, r2, skipidx, axes)

    # Traverse the rest of the kdtree to finish the neighbor list
    find_ball(c_tree.root, nblist, tree_positions, pos, r2, leafnode.leafid,
              skipidx, axes)

    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int find_ball(Node* node,
                   NeighborList nblist,
                   np.float64_t[:, ::1] tree_positions,
                   np.float64_t* pos,
                   np.float64_t r2,
                   uint32_t skipleaf,
                   uint64_t skipidx,
                   axes_range * axes,
                   ) nogil except -1:
    """Traverse the k-d tree to process leaf nodes."""
    if not node.is_leaf:
        if not cull_node_ball(node.less, pos, r2, skipleaf, axes):
            find_ball(node.less, nblist, tree_positions, pos, r2, skipleaf,
                      skipidx, axes)
        if not cull_node_ball(node.greater, pos, r2, skipleaf, axes):
            find_ball(node.greater, nblist, tree_positions, pos, r2, skipleaf,
                      skipidx, axes)
    else:
        if not cull_node_ball(node, pos, r2, skipleaf, axes):
            process_node_points_ball(node, nblist, tree_positions, pos, r2,
                                     skipidx, axes)
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int cull_node_ball(Node* node,
                               np.float64_t* pos,
                               np.float64_t r2,
                               uint32_t skipleaf,
                               axes_range * axes,
                               ) nogil except -1:
    """Check if the node does not intersect with the ball at all."""
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

    return ndist > r2

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int process_node_points_ball(Node* node,
                                         NeighborList nblist,
                                         np.float64_t[:, ::1] positions,
                                         np.float64_t* pos,
                                         np.float64_t r2,
                                         int skipidx,
                                         axes_range * axes,
                                         ) nogil except -1:
    """Add points from the leaf node within the ball to the neighbor list."""
    cdef uint64_t i, k, n
    cdef np.float64_t tpos, sq_dist
    for i in range(node.left_idx, node.left_idx + node.children):
        if i == skipidx:
            continue

        sq_dist = 0.0

        k = axes.start
        while k < axes.stop:
            tpos = positions[i, k] - pos[k]
            sq_dist += tpos*tpos
            k += axes.step

        if (sq_dist < r2):
            nblist.add_pid(sq_dist, i)

    return 0
