"""
This is a recursive function to return a depth-first octree



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython

cdef class position:
    cdef public int output_pos, refined_pos
    def __cinit__(self):
        self.output_pos = 0
        self.refined_pos = 0

cdef class OctreeGrid:
    cdef public object child_indices, fields, left_edges, dimensions, dx
    cdef public int level, offset
    def __cinit__(self,
                  np.ndarray[np.int32_t, ndim=3] child_indices,
                  np.ndarray[np.float64_t, ndim=4] fields,
                  np.ndarray[np.float64_t, ndim=1] left_edges,
                  np.ndarray[np.int32_t, ndim=1] dimensions,
                  np.ndarray[np.float64_t, ndim=1] dx,
                  int level, int offset):
        self.child_indices = child_indices
        self.fields = fields
        self.left_edges = left_edges
        self.dimensions = dimensions
        self.dx = dx
        self.level = level
        self.offset = offset

cdef class OctreeGridList:
    cdef public object grids
    def __cinit__(self, grids):
        self.grids = grids

    def __getitem__(self, int item):
        return self.grids[item]

@cython.boundscheck(False)
def RecurseOctreeDepthFirst(int i_i, int j_i, int k_i,
                            int i_f, int j_f, int k_f,
                            position curpos, int gi, 
                            np.ndarray[np.float64_t, ndim=2] output,
                            np.ndarray[np.int32_t, ndim=1] refined,
                            OctreeGridList grids):
    #cdef int s = curpos
    cdef int i, i_off, j, j_off, k, k_off, ci, fi
    cdef int child_i, child_j, child_k
    cdef OctreeGrid child_grid
    cdef OctreeGrid grid = grids[gi]
    cdef np.ndarray[np.float64_t, ndim=4] fields = grid.fields
    cdef np.ndarray[np.float64_t, ndim=1] leftedges = grid.left_edges
    cdef np.float64_t dx = grid.dx[0]
    cdef np.float64_t child_dx
    cdef np.ndarray[np.float64_t, ndim=1] child_leftedges
    cdef np.float64_t cx, cy, cz
    #here we go over the 8 octants
    #in general however, a mesh cell on this level
    #may have more than 8 children on the next level
    #so we find the int float center (cxyz) of each child cell
    # and from that find the child cell indices
    for i_off in range(i_f):
        i = i_off + i_i #index
        cx = (leftedges[0] + i*dx)
        for j_off in range(j_f):
            j = j_off + j_i
            cy = (leftedges[1] + j*dx)
            for k_off in range(k_f):
                k = k_off + k_i
                cz = (leftedges[2] + k*dx)
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
                    child_grid = grids[ci-grid.offset]
                    child_dx = child_grid.dx[0]
                    child_leftedges = child_grid.left_edges
                    child_i = int((cx - child_leftedges[0])/child_dx)
                    child_j = int((cy - child_leftedges[1])/child_dx)
                    child_k = int((cz - child_leftedges[2])/child_dx)
                    # s = Recurs.....
                    RecurseOctreeDepthFirst(child_i, child_j, child_k, 2, 2, 2,
                                        curpos, ci - grid.offset, output, refined, grids)

@cython.boundscheck(False)
def RecurseOctreeByLevels(int i_i, int j_i, int k_i,
                          int i_f, int j_f, int k_f,
                          np.ndarray[np.int32_t, ndim=1] curpos,
                          int gi, 
                          np.ndarray[np.float64_t, ndim=2] output,
                          np.ndarray[np.int32_t, ndim=2] genealogy,
                          np.ndarray[np.float64_t, ndim=2] corners,
                          OctreeGridList grids):
    cdef np.int32_t i, i_off, j, j_off, k, k_off, ci, fi
    cdef int child_i, child_j, child_k
    cdef OctreeGrid child_grid
    cdef OctreeGrid grid = grids[gi-1]
    cdef int level = grid.level
    cdef np.ndarray[np.int32_t, ndim=3] child_indices = grid.child_indices
    cdef np.ndarray[np.float64_t, ndim=4] fields = grid.fields
    cdef np.ndarray[np.float64_t, ndim=1] leftedges = grid.left_edges
    cdef np.float64_t dx = grid.dx[0]
    cdef np.float64_t child_dx
    cdef np.ndarray[np.float64_t, ndim=1] child_leftedges
    cdef np.float64_t cx, cy, cz
    cdef int cp
    s = None
    for i_off in range(i_f):
        i = i_off + i_i
        cx = (leftedges[0] + i*dx)
        for j_off in range(j_f):
            j = j_off + j_i
            cy = (leftedges[1] + j*dx)
            for k_off in range(k_f):
                k = k_off + k_i
                cz = (leftedges[2] + k*dx)
                cp = curpos[level]
                corners[cp, 0] = cx 
                corners[cp, 1] = cy 
                corners[cp, 2] = cz
                genealogy[curpos[level], 2] = level
                # always output data
                for fi in range(fields.shape[0]):
                    output[cp,fi] = fields[fi,i,j,k]
                ci = child_indices[i,j,k]
                if ci > -1:
                    child_grid = grids[ci-1]
                    child_dx = child_grid.dx[0]
                    child_leftedges = child_grid.left_edges
                    child_i = int((cx-child_leftedges[0])/child_dx)
                    child_j = int((cy-child_leftedges[1])/child_dx)
                    child_k = int((cz-child_leftedges[2])/child_dx)
                    # set current child id to id of next cell to examine
                    genealogy[cp, 0] = curpos[level+1] 
                    # set next parent id to id of current cell
                    genealogy[curpos[level+1]:curpos[level+1]+8, 1] = cp
                    s = RecurseOctreeByLevels(child_i, child_j, child_k, 2, 2, 2,
                                              curpos, ci, output, genealogy,
                                              corners, grids)
                curpos[level] += 1
    return s

