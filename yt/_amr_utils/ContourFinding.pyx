"""
A two-pass contour finding algorithm

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

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

import numpy as np
cimport numpy as np
cimport cython

cdef inline np.int64_t i64max(np.int64_t i0, np.int64_t i1):
    if i0 > i1: return i0
    return i1

cdef inline np.int64_t i64min(np.int64_t i0, np.int64_t i1):
    if i0 < i1: return i0
    return i1

@cython.boundscheck(False)
def construct_boundary_relationships(
        np.ndarray[dtype=np.int64_t, ndim=3] contour_ids):
    # We only look at the boundary and one cell in
    cdef int i, j, nx, ny, nz, offset_i, offset_j, oi, oj
    cdef np.int64_t c1, c2
    tree = []
    nx = contour_ids.shape[0]
    ny = contour_ids.shape[1]
    nz = contour_ids.shape[2]
    # First x-pass
    for i in range(ny):
        for j in range(nz):
            for offset_i in range(3):
                oi = offset_i - 1
                if i == 0 and oi == -1: continue
                if i == ny - 1 and oj == 1: continue
                for offset_j in range(3):
                    oj = offset_j - 1
                    if j == 0 and oj == -1: continue
                    if j == nz - 1 and oj == 1: continue
                    c1 = contour_ids[0, i, j]
                    c2 = contour_ids[1, i + oi, j + oj]
                    if c1 > -1 and c2 > -1:
                        tree.append((i64max(c1,c2), i64min(c1,c2)))
                    c1 = contour_ids[nx-1, i, j]
                    c2 = contour_ids[nx-2, i + oi, j + oj]
                    if c1 > -1 and c2 > -1:
                        tree.append((i64max(c1,c2), i64min(c1,c2)))
    # Now y-pass
    for i in range(nx):
        for j in range(nz):
            for offset_i in range(3):
                oi = offset_i - 1
                if i == 0 and oi == -1: continue
                if i == nx - 1 and oj == 1: continue
                for offset_j in range(3):
                    oj = offset_j - 1
                    if j == 0 and oj == -1: continue
                    if j == nz - 1 and oj == 1: continue
                    c1 = contour_ids[i, 0, j]
                    c2 = contour_ids[i + oi, 1, j + oj]
                    if c1 > -1 and c2 > -1:
                        tree.append((i64max(c1,c2), i64min(c1,c2)))
                    c1 = contour_ids[i, ny-1, j]
                    c2 = contour_ids[i + oi, ny-2, j + oj]
                    if c1 > -1 and c2 > -1:
                        tree.append((i64max(c1,c2), i64min(c1,c2)))
    for i in range(nx):
        for j in range(ny):
            for offset_i in range(3):
                oi = offset_i - 1
                if i == 0 and oi == -1: continue
                if i == nx - 1 and oj == 1: continue
                for offset_j in range(3):
                    oj = offset_j - 1
                    if j == 0 and oj == -1: continue
                    if j == ny - 1 and oj == 1: continue
                    c1 = contour_ids[i, j, 0]
                    c2 = contour_ids[i + oi, j + oj, 1]
                    if c1 > -1 and c2 > -1:
                        tree.append((i64max(c1,c2), i64min(c1,c2)))
                    c1 = contour_ids[i, j, nz-1]
                    c2 = contour_ids[i + oi, j + oj, nz-2]
                    if c1 > -1 and c2 > -1:
                        tree.append((i64max(c1,c2), i64min(c1,c2)))
    return tree
