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
from cykdtree.kdtree cimport PyKDTree, KDTree, Node, uint64_t
from cython.operator cimport dereference as deref

cdef extern from "<algorithm>" namespace "std" nogil:
    void sort[Iter](Iter first, Iter last)

from libcpp.vector cimport vector
from libc.math cimport sqrt

from yt.funcs import get_pbar

@cython.boundscheck(False)
@cython.wraparound(False)
def generate_smoothing_length(np.float64_t[:, ::1] input_positions,
                              PyKDTree kdtree,
                              int n_neighbors):
    """Calculate array of distances to the nth nearest neighbor

    Parameters
    ----------

    input_positions: arrays of floats with shape (nparticles, 3)
        The positions of particles. Current assumed to be 3D postions.
    kdtree: A PyKDTree instance
        A kdtree to do nearest neighbors searches with
    n_neighbors: The neighbor number to calculate the distance to

    Returns
    -------

    smoothing_lengths: arrays of flots with shape (nparticles, )
        The calculated smoothing lengths

    """
    cdef KDTree* c_tree = kdtree._tree
    cdef Node* leafnode
    cdef Node* neighbor
    cdef uint64_t n_nearby
    cdef vector[np.float64_t] squared_distances
    cdef vector[uint64_t] nearby_indices
    cdef int n_particles = input_positions.shape[0]
    cdef np.float64_t[:] smoothing_length = np.empty(n_particles)
    cdef np.float64_t tpos
    cdef np.float64_t* nearby_position
    cdef np.float64_t* pos
    cdef np.float64_t* positions = &input_positions[0, 0]
    cdef uint64_t neighbor_id
    cdef int i, j, k
    pbar = get_pbar("Calculating smoothing length", n_particles)
    with nogil:
        for i in range(n_particles):
            if i % 1000 == 0:
                with gil:
                    pbar.update(i)
                    PyErr_CheckSignals()
            # This section is optimized. Here we're using pointer arithmetic
            # to get a reference to the place in input_positions where the
            # position of the particle under consideration is stored
            pos = positions + 3*i

            # Search the tree for the node containing the position under
            # consideration. Fine neighbor nodes and determine the total
            # number of particles in all the nodes under consideration
            leafnode = c_tree.search(&pos[0])
            n_nearby = leafnode.children
            for neighbor_id in leafnode.all_neighbors:
                n_nearby = n_nearby + c_tree.leaves[neighbor_id].children

            # Find indices into the particle position array for the list of
            # potential nearest neighbors
            nearby_indices = vector[uint64_t](n_nearby)
            for j in range(leafnode.children):
                nearby_indices[j] = c_tree.all_idx[leafnode.left_idx+j]
            for neighbor_id in leafnode.all_neighbors:
                neighbor = c_tree.leaves[neighbor_id]
                for k in range(neighbor.children):
                    nearby_indices[j+k] = c_tree.all_idx[neighbor.left_idx+k]
                j += neighbor.children

            # Calculate the squared distances to all of the particles in
            # the neighbor list
            squared_distances = vector[np.float64_t](n_nearby)
            for j in range(n_nearby):
                nearby_position = positions + 3*nearby_indices[j]
                for k in range(3):
                    tpos = nearby_position[k] - pos[k]
                    squared_distances[j] += tpos*tpos

            # Sort the squared distances and find the nth entry, this is the
            # nth nearest neighbor for particle i. Take the square root of
            # the squared distance to the nth neighbor to find the smoothing
            # length.
            sort[vector[np.float64_t].iterator](
                squared_distances.begin(), squared_distances.end())
            smoothing_length[i] = sqrt(squared_distances[n_neighbors])
    pbar.finish()
    return np.array(smoothing_length)
