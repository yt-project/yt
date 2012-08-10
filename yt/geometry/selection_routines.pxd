"""
Geometry selection routine imports.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

cdef struct Oct

cdef class SelectorObject:
    cdef void recursively_select_octs(self, Oct *root,
                        np.float64_t pos[3], np.float64_t dds[3],
                        np.ndarray[np.uint8_t, ndim=2] mask,
                        int level = ?)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3],
                         int eterm[3]) nogil
    cdef void set_bounds(self,
                         np.float64_t left_edge[3], np.float64_t right_edge[3],
                         np.float64_t dds[3], int ind[3][2], int *check)
