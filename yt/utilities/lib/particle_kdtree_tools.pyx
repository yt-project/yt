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
from cython.operator cimport dereference as deref

cdef extern from "<algorithm>" namespace "std" nogil:
    void sort[Iter](Iter first, Iter last)

from libcpp.vector cimport vector
from libc.math cimport sqrt

from yt.funcs import get_pbar

cdef int CHUNKSIZE = 4096

@cython.boundscheck(True)
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
    cdef Node* node
    cdef vector[np.float64_t] squared_distances
    cdef vector[uint64_t] nearby_indices
    cdef vector[uint32_t] nearby_ids
    cdef int n_particles = input_positions.shape[0]
    cdef np.float64_t[:] smoothing_length = np.empty(n_particles)
    cdef np.float64_t tpos, furthest_distance, sq_dist
    cdef uint64_t neighbor_id
    cdef int i, j, k, l, n_kept, skip
    pbar = get_pbar("Generate smoothing length", n_particles)
    for i in range(n_particles):
        if i % CHUNKSIZE == 0:
            pbar.update(i-1)
            PyErr_CheckSignals()

        # Search the tree for the node containing the position under
        # consideration. Fine neighbor nodes and determine the total
        # number of particles in all the nodes under consideration
        leafnode = c_tree.search(&input_positions[i, 0])

        # Find indices into the particle position array for the list of
        # potential nearest neighbors
        nearby_ids = leafnode.all_neighbors
        nearby_ids.push_back(leafnode.leafid)

        squared_distances = vector[np.float64_t]()
        furthest_distance = 0
        n_kept = 0
        for j in range(nearby_ids.size() - 1, -1, -1):
            node = c_tree.leaves[nearby_ids[j]]

            nearby_indices = vector[uint64_t]()
            for k in range(node.children):
                nearby_indices.push_back(c_tree.all_idx[node.left_idx + k])

            for l in range(nearby_indices.size()):
                skip = 0
                sq_dist = 0
                for k in range(3):
                    tpos = (input_positions[nearby_indices[l], k] -
                            input_positions[i, k])
                    sq_dist += tpos*tpos
                    if (n_kept > n_neighbors and sq_dist > furthest_distance):
                        skip = 1
                        break
                if skip:
                    continue
                squared_distances.push_back(sq_dist)
                n_kept += 1
                if squared_distances[j] > furthest_distance:
                    furthest_distance = squared_distances[j]

        # Sort the squared distances and find the nth entry, this is the
        # nth nearest neighbor for particle i. Take the square root of
        # the squared distance to the nth neighbor to find the smoothing
        # length.
        sort[vector[np.float64_t].iterator](
            squared_distances.begin(), squared_distances.end())
        smoothing_length[i] = sqrt(squared_distances[n_neighbors])
    pbar.update(n_particles-1)
    pbar.finish()
    return np.array(smoothing_length)
