"""
Make a fake octree, deposit particle at every leaf




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from libc.stdlib cimport malloc, free, rand, RAND_MAX
cimport numpy as np
from oct_visitors cimport cind
import numpy as np
cimport cython

from oct_container cimport Oct, RAMSESOctreeContainer

# Create a balanced octree by a random walk that recursively
# subdivides
def create_fake_octree(RAMSESOctreeContainer oct_handler,
                       long max_noct,
                       long max_level,
                       np.ndarray[np.int32_t, ndim=1] ndd,
                       np.ndarray[np.float64_t, ndim=1] dle,
                       np.ndarray[np.float64_t, ndim=1] dre,
                       float fsubdivide):
    cdef int[3] dd #hold the octant index
    cdef int[3] ind #hold the octant index
    cdef long i
    cdef long cur_leaf = 0
    cdef np.ndarray[np.uint8_t, ndim=2] mask
    for i in range(3):
        ind[i] = 0
        dd[i] = ndd[i]
    oct_handler.allocate_domains([max_noct])
    parent = oct_handler.next_root(1, ind)
    parent.domain = 1
    cur_leaf = 8 #we've added one parent...
    mask = np.ones((max_noct,8),dtype='uint8')
    while oct_handler.domains[0].n_assigned < max_noct:
        print "root: nocts ", oct_handler.domains[0].n_assigned
        cur_leaf = subdivide(oct_handler, parent, ind, dd, cur_leaf, 0,
                             max_noct, max_level, fsubdivide, mask)
    return cur_leaf


cdef long subdivide(RAMSESOctreeContainer oct_handler,
                    Oct *parent,
                    int ind[3], int dd[3],
                    long cur_leaf, long cur_level,
                    long max_noct, long max_level, float fsubdivide,
                    np.ndarray[np.uint8_t, ndim=2] mask):
    print "child", parent.file_ind, ind[0], ind[1], ind[2], cur_leaf, cur_level
    cdef int ddr[3]
    cdef int ii
    cdef long i
    cdef float rf #random float from 0-1
    if cur_level >= max_level:
        return cur_leaf
    if oct_handler.domains[0].n_assigned >= max_noct:
        return cur_leaf
    for i in range(3):
        ind[i] = <int> ((rand() * 1.0 / RAND_MAX) * dd[i])
        ddr[i] = 2
    rf = rand() * 1.0 / RAND_MAX
    if rf > fsubdivide:
        ii = cind(ind[0], ind[1], ind[2])
        if parent.children[ii] == NULL:
            cur_leaf += 7
        oct = oct_handler.next_child(1, ind, parent)
        oct.domain = 1
        cur_leaf = subdivide(oct_handler, oct, ind, ddr, cur_leaf,
                             cur_level + 1, max_noct, max_level,
                             fsubdivide, mask)
    return cur_leaf
