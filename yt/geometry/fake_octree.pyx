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

# Defined only by N leaves
# Randomly decide if a branch should be subdivide, recurse one level if so
# Once done, create a position array of len(leafes) with smoothing lengths = oct_size

# Note that with this algorithm the octree won't be balanced once you hit
# the maximum number of desired leaves

# Use next_child(domain, int[3] octant, Oct parent)

def create_fake_octree(RAMSESOctreeContainer oct_handler,
                       long noct,
                       long max_level,
                       np.ndarray[np.int32_t, ndim=1] ndd,
                       np.ndarray[np.float64_t, ndim=1] dle,
                       np.ndarray[np.float64_t, ndim=1] dre,
                       float fsubdivide):
    cdef int[3] ind #hold the octant index
    cdef int[3] dd #hold the octant index
    cdef long i
    cdef long cur_noct = 0
    for i in range(3):
        ind[i] = 0
        dd[i] = ndd[i]
    assert dd[0]*dd[1]*dd[2] <= noct
    print 'starting'
    print ind[0], ind[1], ind[2]
    print 'allocate'
    print noct
    oct_handler.allocate_domains([noct])
    print 'n_assigned', oct_handler.domains[0].n_assigned
    print 'parent'
    parent = oct_handler.next_root(oct_handler.max_domain, ind)
    print 'subdiv'
    while oct_handler.domains[0].n_assigned < noct:
        cur_noct = subdivide(oct_handler,ind, dd, parent, 0, 0, noct,
                  max_level, fsubdivide)

cdef long subdivide(RAMSESOctreeContainer oct_handler, int ind[3], 
               int dd[3],
               Oct *parent, long cur_level, long cur_noct,
               long noct, long max_level, float fsubdivide):
    print cur_level, ' n_assigned ', oct_handler.domains[0].n_assigned, 
    print ' n', oct_handler.domains[0].n
    cdef int ddr[3]
    cdef long i,j,k
    cdef float rf #random float from 0-1
    if cur_level >= max_level: 
        return cur_noct
    if oct_handler.domains[0].n_assigned >= noct: 
        return cur_noct
    for i in range(3):
        ind[i] = <int> ((rand() * 1.0 / RAND_MAX) * dd[i])
        ddr[i] = 2
    rf = rand() * 1.0 / RAND_MAX
    if rf > fsubdivide:
        #this will mark the octant ind as subdivided
        oct = oct_handler.next_child(1, ind, parent)
        subdivide(oct_handler, ind, ddr, oct, cur_level + 1, 
                  cur_noct+ 1, noct, max_level, fsubdivide)
    return cur_noct
