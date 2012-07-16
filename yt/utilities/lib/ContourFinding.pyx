"""
A two-pass contour finding algorithm

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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
from libc.stdlib cimport malloc, free

cdef extern from "math.h":
    double fabs(double x)

cdef inline np.int64_t i64max(np.int64_t i0, np.int64_t i1):
    if i0 > i1: return i0
    return i1

cdef inline np.int64_t i64min(np.int64_t i0, np.int64_t i1):
    if i0 < i1: return i0
    return i1

cdef extern from "union_find.h":
    ctypedef struct forest_node:
        void *value
        forest_node *parent
        int rank

    forest_node* MakeSet(void* value)
    void Union(forest_node* node1, forest_node* node2)
    forest_node* Find(forest_node* node)

ctypedef struct CellIdentifier:
    np.int64_t hindex
    int level

cdef class GridContourContainer:
    cdef np.int64_t dims[3]
    cdef np.int64_t start_indices[3]
    cdef forest_node **join_tree
    cdef np.int64_t ncells

    def __init__(self, dimensions, indices):
        cdef int i
        self.ncells = 1
        for i in range(3):
            self.ncells *= dimensions[i]
            self.dims[i] = dimensions[i]
            self.start_indices[i] = indices[i]
        self.join_tree = <forest_node **> malloc(sizeof(forest_node) 
                                                 * self.ncells)
        for i in range(self.ncells): self.join_tree[i] = NULL

    def __dealloc__(self):
        cdef int i
        for i in range(self.ncells):
            if self.join_tree[i] != NULL: free(self.join_tree[i])
        free(self.join_tree)

    #def construct_join_tree(self,
    #            np.ndarray[np.float64_t, ndim=3] field,
    #            np.ndarray[np.bool_t, ndim=3] mask):
    #    # This only looks at the components of the grid that are actually
    #    # inside this grid -- boundary conditions are handled later.
    #    pass

#@cython.boundscheck(False)
#@cython.wraparound(False)
def construct_boundary_relationships(
        np.ndarray[dtype=np.int64_t, ndim=3] contour_ids):
    # We only look at the boundary and one cell in
    cdef int i, j, nx, ny, nz, offset_i, offset_j, oi, oj
    cdef np.int64_t c1, c2
    nx = contour_ids.shape[0]
    ny = contour_ids.shape[1]
    nz = contour_ids.shape[2]
    # We allocate an array of fixed (maximum) size
    cdef int s = (ny*nx + nx*nz + ny*nz - 2) * 18
    cdef np.ndarray[np.int64_t, ndim=2] tree = np.zeros((s, 2), dtype="int64")
    cdef int ti = 0
    # First x-pass
    for i in range(ny):
        for j in range(nz):
            for offset_i in range(3):
                oi = offset_i - 1
                if i == 0 and oi == -1: continue
                if i == ny - 1 and oi == 1: continue
                for offset_j in range(3):
                    oj = offset_j - 1
                    if j == 0 and oj == -1: continue
                    if j == nz - 1 and oj == 1: continue
                    c1 = contour_ids[0, i, j]
                    c2 = contour_ids[1, i + oi, j + oj]
                    if c1 > -1 and c2 > -1:
                        tree[ti,0] = i64max(c1,c2)
                        tree[ti,1] = i64min(c1,c2)
                        ti += 1
                    c1 = contour_ids[nx-1, i, j]
                    c2 = contour_ids[nx-2, i + oi, j + oj]
                    if c1 > -1 and c2 > -1:
                        tree[ti,0] = i64max(c1,c2)
                        tree[ti,1] = i64min(c1,c2)
                        ti += 1
    # Now y-pass
    for i in range(nx):
        for j in range(nz):
            for offset_i in range(3):
                oi = offset_i - 1
                if i == 0 and oi == -1: continue
                if i == nx - 1 and oi == 1: continue
                for offset_j in range(3):
                    oj = offset_j - 1
                    if j == 0 and oj == -1: continue
                    if j == nz - 1 and oj == 1: continue
                    c1 = contour_ids[i, 0, j]
                    c2 = contour_ids[i + oi, 1, j + oj]
                    if c1 > -1 and c2 > -1:
                        tree[ti,0] = i64max(c1,c2)
                        tree[ti,1] = i64min(c1,c2)
                        ti += 1
                    c1 = contour_ids[i, ny-1, j]
                    c2 = contour_ids[i + oi, ny-2, j + oj]
                    if c1 > -1 and c2 > -1:
                        tree[ti,0] = i64max(c1,c2)
                        tree[ti,1] = i64min(c1,c2)
                        ti += 1
    for i in range(nx):
        for j in range(ny):
            for offset_i in range(3):
                oi = offset_i - 1
                if i == 0 and oi == -1: continue
                if i == nx - 1 and oi == 1: continue
                for offset_j in range(3):
                    oj = offset_j - 1
                    if j == 0 and oj == -1: continue
                    if j == ny - 1 and oj == 1: continue
                    c1 = contour_ids[i, j, 0]
                    c2 = contour_ids[i + oi, j + oj, 1]
                    if c1 > -1 and c2 > -1:
                        tree[ti,0] = i64max(c1,c2)
                        tree[ti,1] = i64min(c1,c2)
                        ti += 1
                    c1 = contour_ids[i, j, nz-1]
                    c2 = contour_ids[i + oi, j + oj, nz-2]
                    if c1 > -1 and c2 > -1:
                        tree[ti,0] = i64max(c1,c2)
                        tree[ti,1] = i64min(c1,c2)
                        ti += 1
    return tree[:ti,:]

cdef inline int are_neighbors(
            np.float64_t x1, np.float64_t y1, np.float64_t z1,
            np.float64_t dx1, np.float64_t dy1, np.float64_t dz1,
            np.float64_t x2, np.float64_t y2, np.float64_t z2,
            np.float64_t dx2, np.float64_t dy2, np.float64_t dz2,
        ):
    # We assume an epsilon of 1e-15
    if fabs(x1-x2) > 0.5*(dx1+dx2): return 0
    if fabs(y1-y2) > 0.5*(dy1+dy2): return 0
    if fabs(z1-z2) > 0.5*(dz1+dz2): return 0
    return 1

@cython.boundscheck(False)
@cython.wraparound(False)
def identify_field_neighbors(
            np.ndarray[dtype=np.float64_t, ndim=1] field,
            np.ndarray[dtype=np.float64_t, ndim=1] x,
            np.ndarray[dtype=np.float64_t, ndim=1] y,
            np.ndarray[dtype=np.float64_t, ndim=1] z,
            np.ndarray[dtype=np.float64_t, ndim=1] dx,
            np.ndarray[dtype=np.float64_t, ndim=1] dy,
            np.ndarray[dtype=np.float64_t, ndim=1] dz,
        ):
    # We assume this field is pre-jittered; it has no identical values.
    cdef int outer, inner, N, added
    cdef np.float64_t x1, y1, z1, dx1, dy1, dz1
    N = field.shape[0]
    #cdef np.ndarray[dtype=np.object_t] joins
    joins = [[] for outer in range(N)]
    #joins = np.empty(N, dtype='object')
    for outer in range(N):
        if (outer % 10000) == 0: print outer, N
        x1 = x[outer]
        y1 = y[outer]
        z1 = z[outer]
        dx1 = dx[outer]
        dy1 = dy[outer]
        dz1 = dz[outer]
        this_joins = joins[outer]
        added = 0
        # Go in reverse order
        for inner in range(outer, 0, -1):
            if not are_neighbors(x1, y1, z1, dx1, dy1, dz1,
                                 x[inner], y[inner], z[inner],
                                 dx[inner], dy[inner], dz[inner]):
                continue
            # Hot dog, we have a weiner!
            this_joins.append(inner)
            added += 1
            if added == 26: break
    return joins

@cython.boundscheck(False)
@cython.wraparound(False)
def extract_identified_contours(int max_ind, joins):
    cdef int i
    contours = []
    for i in range(max_ind + 1): # +1 to get to the max_ind itself
        contours.append(set([i]))
        if len(joins[i]) == 0:
            continue
        proto_contour = [i]
        for j in joins[i]:
            proto_contour += contours[j]
        proto_contour = set(proto_contour)
        for j in proto_contour:
            contours[j] = proto_contour
    return contours

@cython.boundscheck(False)
@cython.wraparound(False)
def update_joins(joins, np.ndarray[np.int64_t, ndim=1] contour_ids):
    cdef np.int64_t new, old, i, oi
    cdef int n, on
    cdef np.ndarray[np.int64_t, ndim=1] old_set
    #print contour_ids.shape[0]
    n = contour_ids.shape[0]
    for new, old_set in joins:
        #print new
        on = old_set.shape[0]
        for i in range(n):
            for oi in range(on):
                old = old_set[oi]
                if contour_ids[i] == old: contour_ids[i] = new
