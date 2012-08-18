"""
Oct definitions file

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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

cimport numpy as np

cdef struct ParticleArrays

cdef struct Oct
cdef struct Oct:
    np.int64_t ind          # index
    np.int64_t local_ind
    np.int64_t domain       # (opt) addl int index
    np.int64_t pos[3]       # position in ints
    np.int8_t level
    ParticleArrays *sd
    Oct *children[2][2][2]
    Oct *parent

cdef struct OctAllocationContainer
cdef struct OctAllocationContainer:
    np.int64_t n
    np.int64_t n_assigned
    np.int64_t offset
    OctAllocationContainer *next
    Oct *my_octs

cdef class OctreeContainer:
    cdef OctAllocationContainer *cont
    cdef Oct ****root_mesh
    cdef int nn[3]
    cdef np.float64_t DLE[3], DRE[3]
    cdef public int nocts
    cdef public int max_domain

cdef class RAMSESOctreeContainer(OctreeContainer):
    cdef OctAllocationContainer **domains

cdef struct ParticleArrays:
    np.float64_t **pos
    np.int64_t *domain_id
    np.int64_t np
