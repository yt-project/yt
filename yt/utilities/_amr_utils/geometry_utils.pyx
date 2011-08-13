"""
Simple integrators for the radiative transfer equation

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

import numpy as np
cimport numpy as np
cimport cython
from stdlib cimport malloc, free, abs

def ortho_ray_grids(dobj, np.ndarray[np.float64_t, ndim=2] left_edges,
                          np.ndarray[np.float64_t, ndim=2] right_edges):
    cdef int i
    cdef int ng = left_edges.shape[0]
    cdef int px_ax = dobj.px_ax
    cdef int py_ax = dobj.py_ax
    cdef np.float64_t px = dobj.px
    cdef np.float64_t py = dobj.py
    cdef np.ndarray[np.int32_t, ndim=1] gridi = np.zeros(ng, dtype='int32_t')
    for i in range(ng):
        if (    (px >= left_edges[i, px])
            and (px < right_edges[i, px])
            and (py >= left_edges[i, py])
            and (py < right_edges[i, py])):
            gridi[i] = 1
    return gridi

def ray_grids(dobj, np.ndarray[np.float64_t, ndim=2] left_edges,
                    np.ndarray[np.float64_t, ndim=2] right_edges):
    cdef int i, ax
    cdef int i1, i2
    cdef int ng = left_edges.shape[0]
    cdef np.ndarray[np.int32_t, ndim=1] gridi = np.zeros(ng, dtype='int32_t')
    cdef np.float64_t vs[3], t, p0[3], p1[3], v[3]
    for i in range(3):
        p0[i] = dobj.start_point[i]
        p1[i] = dobj.end_point[i]
        v[i] = dobj.vec[i]
    # We check first to see if at any point, the ray intersects a grid face
    for gi in range(ng):
        for ax in range(3):
            i1 = (ax+1) % 3
            i2 = (ax+2) % 3
            t = (left_edges[gi,ax] - p0[ax])/v[ax]
            for i in range(3)
                vs[i] = t * v[i] + p0[i]
            if left_edges[gi,i1] <= vs[i1] and \
               right_edges[gi,i1] >= vs[i1] and \
               left_edges[gi,i2] <= vs[i2] and \
               right_edges[gi,i2] >= vs[i2]:
                gridi[gi] = 1
                break
            t = (right_edges[gi,ax] - p0[ax])/v[ax]
            for i in range(3)
                vs[i] = t * v[i] + p0[i]
            if left_edges[gi,i1] <= vs[i1] and \
               right_edges[gi,i1] >= vs[i1] and \
               left_edges[gi,i2] <= vs[i2] and \
               right_edges[gi,i2] >= vs[i2]:
                gridi[gi] = 1
                break
        if gridi[gi] = 1: continue
        # if the point is fully enclosed, we count the grid
        if left_edges[gi,0] <= p0[0] and \
           right_edges[gi,0] >= p0[0]:
           left_edges[gi,1] <= p0[1] and \
           right_edges[gi,1] >= p0[1]:
           left_edges[gi,2] <= p0[2] and \
           right_edges[gi,2] >= p0[2]:
            gridi[gi] = 1
            continue
        if left_edges[gi,0] <= p1[0] and \
           right_edges[gi,0] >= p1[0]:
           left_edges[gi,1] <= p1[1] and \
           right_edges[gi,1] >= p1[1]:
           left_edges[gi,2] <= p1[2] and \
           right_edges[gi,2] >= p1[2]:
            gridi[gi] = 1
            continue
    return gridi
