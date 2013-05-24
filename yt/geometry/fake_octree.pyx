"""
Make a fake octree, deposit particle at every leaf

Author: Christopher Moody <chris.e.moody@gmail.com>
Affiliation: UC Santa Cruz
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2013 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from libc.stdlib cimport malloc, free, rand, RAND_MAX
cimport numpy as np
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
    cdef long i,j,k
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
        if parent.children[ind[0]][ind[1]][ind[2]] == NULL:
            cur_leaf += 7 
        oct = oct_handler.next_child(1, ind, parent)
        oct.domain = 1
        cur_leaf = subdivide(oct_handler, oct, ind, ddr, cur_leaf, 
                             cur_level + 1, max_noct, max_level, 
                             fsubdivide, mask)
    return cur_leaf
