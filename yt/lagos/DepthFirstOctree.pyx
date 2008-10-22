"""
This is a recursive function to return a depth-first octree

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

cdef class position:
    cdef public int output_pos, refined_pos

cdef class OctreeGrid:
    cdef public object child_indices, fields, left_edges, dimensions, dx
    def __cinit__(self,
                  np.ndarray[np.int32_t, ndim=3] child_indices,
                  np.ndarray[np.float64_t, ndim=4] fields,
                  np.ndarray[np.float64_t, ndim=1] left_edges,
                  np.ndarray[np.int32_t, ndim=1] dimensions,
                  np.ndarray[np.float64_t, ndim=1] dx):
        self.child_indices = child_indices
        self.fields = fields
        self.left_edges = left_edges
        self.dimensions = dimensions
        self.dx = dx

cdef class OctreeGridList:
    cdef public object grids
    def __cinit__(self, grids):
        self.grids = grids

    def __getitem__(self, int item):
        return self.grids[item]

@cython.boundscheck(False)
def WalkRootgrid(np.ndarray[np.float64_t, ndim=2] output,
                 np.ndarray[np.int32_t, ndim=1] refined,
                 OctreeGridList grids, int pi, int s = 0, int r = 0):
    """
    This function only gets called on a 'root grid' -- a base grid
    of the simulation we are converting.  It will call a recursive walker.
    """
    cdef int i, j, k, gi, fi
    cdef int child_i, child_j, child_k
    cdef OctreeGrid child_grid
    cdef OctreeGrid grid = grids[pi-1]
    cdef np.ndarray[np.int32_t, ndim=3] child_indices = grid.child_indices
    cdef np.ndarray[np.int32_t, ndim=1] dimensions = grid.dimensions
    cdef np.ndarray[np.float64_t, ndim=4] fields = grid.fields
    cdef np.ndarray[np.float64_t, ndim=1] leftedges = grid.left_edges
    cdef np.ndarray[np.float64_t, ndim=1] dx = grid.dx
    cdef np.ndarray[np.float64_t, ndim=1] child_dx
    cdef np.ndarray[np.float64_t, ndim=1] child_leftedges
    cdef position curpos
    curpos.output_pos = s
    curpos.refined_pos = r
    for k in range(dimensions[2]):
        print k
        for j in range(dimensions[1]):
            for i in range(dimensions[0]):
                gi = child_indices[i,j,k]
                if gi == -1:
                    for fi in range(fields.shape[0]):
                        output[curpos.output_pos,fi] = fields[fi,i,j,k]
                    refined[curpos.refined_pos] = 0
                    curpos.output_pos += 1
                    curpos.refined_pos += 1
                else:
                    refined[curpos.refined_pos] = 1
                    curpos.refined_pos += 1
                    child_grid = grids[gi-1]
                    child_dx = child_grid.dx
                    child_leftedges = child_grid.left_edges
                    child_i = ((leftedges[0] + i * dx) - child_leftedges[0])/child_dx
                    child_j = ((leftedges[1] + j * dx) - child_leftedges[1])/child_dx
                    child_k = ((leftedges[2] + k * dx) - child_leftedges[2])/child_dx
                    RecurseOctree(child_i, child_j, child_k, curpos, gi, pi, output, refined, grids)

@cython.boundscheck(False)
def RecurseOctree(int i_i, int j_i, int k_i,
                  position curpos, int gi, int pi,
                  np.ndarray[np.float64_t, ndim=2] output,
                  np.ndarray[np.int32_t, ndim=1] refined,
                  OctreeGridList grids):
    cdef int i, i_off, j, j_off, k, k_off, ci, fi
    cdef child_i, child_j, child_k
    cdef OctreeGrid child_grid
    cdef OctreeGrid grid = grids[gi-1]
    cdef np.ndarray[np.int32_t, ndim=3] child_indices = grid.child_indices
    cdef np.ndarray[np.int32_t, ndim=1] dimensions = grid.dimensions
    cdef np.ndarray[np.float64_t, ndim=4] fields = grid.fields
    cdef np.ndarray[np.float64_t, ndim=1] leftedges = grid.left_edges
    cdef np.ndarray[np.float64_t, ndim=1] dx = grid.dx
    cdef np.ndarray[np.float64_t, ndim=1] child_dx
    cdef np.ndarray[np.float64_t, ndim=1] child_leftedges
    # Not sure how to get around this
    for k_off in range(2):
        k = k_off + k_i
        for j_off in range(2):
            j = j_off + j_i
            for i_off in range(2):
                i = i_off + i_i
                ci = grid.child_indices[i,j,k]
                if ci == -1:
                    for fi in range(fields.shape[0]):
                        output[curpos.output_pos,fi] = fields[fi,i,j,k]
                    refined[curpos.refined_pos] = 0
                    curpos.output_pos += 1
                    curpos.refined_pos += 1
                else:
                    refined[curpos.refined_pos] = 1
                    curpos.refined_pos += 1
                    child_grid = grids[ci-1]
                    child_dx = child_grid.dx
                    child_leftedges = child_grid.left_edges
                    child_i = ((leftedges[0] + i * dx) - child_leftedges[0])/child_dx
                    child_j = ((leftedges[1] + j * dx) - child_leftedges[1])/child_dx
                    child_k = ((leftedges[2] + k * dx) - child_leftedges[2])/child_dx
                    s = RecurseOctree(child_i, child_j, child_k, curpos, ci, gi, output, refined, grids)
    return s
