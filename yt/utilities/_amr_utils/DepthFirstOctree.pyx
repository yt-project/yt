"""
This is a recursive function to return a depth-first octree

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
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
    def __cinit__(self):
        self.output_pos = 0
        self.refined_pos = 0

cdef class hilbert_state:
#    """
#    From libpjutil/hilbert.h:
#    This class represents the extra information associated with an
#    octree cell to be able to perform a Hilbert ordering. This
#    information consists of a permutation of the dimensions and the
#    direction of the 3 axes, which can be encoded in one byte. For an
#    octree traversal stack this information needs to be saved at each
#    level, so it's good to have a compact representation. 
#
#    The data is stored as follows: Bits 0-1 stores the first
#    dimension, 2-3 the second. Because we know it is a permutation of
#    012, there is no need to store the third dimension. Then bits 4-6
#    are the signs of the three axes.
#
#    Apparently there is no way to encode a uint8_t literal except as a
#    character constant, hence the use of those.
#    """
# These assignments are from 2.7 of BG200
# vertex > 1st dim 2nd dim 3rd dim
# 1 > +z+x+y 
# 2 > +y+z+x
# 3 > +y+x+z
# 4 > -z-y+x 
# 5 > +z-y-x
# 6 > +y+x+z
# 7 > +y-z-x
# 8 > -z+x-y
    cdef public int dima,dimb,dimc,signa,signb,signc
    #cdef np.ndarray[np.int32,ndim =1] a = na.zeros(3)
    #cdef np.ndarray[np.int32,ndim =1] b = na.zeros(3)
    #cdef np.ndarray[np.int32,ndim =1] c = na.zeros(3)
    #cdef np.ndarray[np.int32,ndim =1] d = na.zeros(3)
    #cdef np.ndarray[np.int32,ndim =1] e = na.zeros(3)
    #cdef np.ndarray[np.int32,ndim =1] f = na.zeros(3)
    #cdef np.ndarray[np.int32,ndim =1] g = na.zeros(3)
    #cdef np.ndarray[np.int32,ndim =1] h = na.zeros(3)

    cdef np.ndarray[np.int64,ndim =1] oct = na.zeros(3)
    cdef np.ndarray[np.int64,ndim =1] noct = na.zeros(3)


    #a[0],a[1],a[2] = 0,0,0
    #b[0],b[1],b[2] = 0,1,0
    #c[0],c[1],c[2] = 1,1,0
    #d[0],d[1],d[2] = 1,0,0
    #e[0],e[1],e[2] = 1,0,1
    #f[0],f[1],f[2] = 1,1,1
    #g[0],g[1],g[2] = 0,1,1
    #h[0],h[1],h[2] = 0,0,1

    def __cinit__(int parent_octant):
        self.swap_by_octant(parent_octant)
        self.signc = signa*signb
        self.dimc = 3-dima-dimb
    
    def swap_by_octant(int octant): 
        if octant==1: 
            self.dima = 2
            self.dimb = 0
            self.signa = 1
            self.signb = 1
        if octant==2: 
            self.dima = 1
            self.dimb = 2
            self.signa = 1
            self.signb = 1
        if octant==3:
            self.dima = 1
            self.dimb = 0
            self.signa = 1
            self.signb = 1
        if octant==4: 
            self.dima = 2
            self.dimb = 1
            self.signa = -1
            self.signb = -1
        if octant==5: 
            self.dima = 2
            self.dimb = 1
            self.signa = 1
            self.signb = -1
        if octant==6: 
            self.dima = 1
            self.dimb = 0
            self.signa = 1
            self.signb = 1
        if octant==7: 
            self.dima = 1
            self.dimb = 2
            self.signa = 1
            self.signb = -1
        if octant==8: 
            self.dima = 2
            self.dimb = 0
            self.signa = -1
            self.signb = 1

    def __iter__(self):
        return self.next_hilbert()

    def next(self):
        #yield the next cell in this oct
        
        #as/descend the first dimension
        # the second dim
        #reverse the first
        #climb the third
        oct[self.dima] = 0 if self.signx>0 else 1
        oct[self.dimb] = 0 if self.signy>0 else 1
        oct[self.dimc] = 0 if self.signz>0 else 1
        yield oct
        oct[self.dima] += self.signa; yield oct
        oct[self.dimb] += self.signb; yield oct
        oct[self.dima] -= self.signa; yield oct
        oct[self.dimc] += self.signc; yield oct
        oct[self.dima] += self.signa; yield oct
        oct[self.dimb] -= self.signb; yield oct
        oct[self.dimb] -= self.signa; return oct

    def next_hilbert(self):
        noct = self.next()
        return noct, hilbert_state(noct)

@cython.boundscheck(False)
def RecurseOctreeDepthFirstHilbert(int i_i, int j_i, int k_i,
                            int i_f, int j_f, int k_f,
                            position curpos, int gi, 
                            hilbert_state hs,
                            np.ndarray[np.float64_t, ndim=3] output,
                            np.ndarray[np.int64_t, ndim=1] refined,
                            OctreeGridList grids):
    #cdef int s = curpos
    cdef int i, i_off, j, j_off, k, k_off, ci, fi
    cdef int child_i, child_j, child_k
    cdef OctreeGrid child_grid
    cdef OctreeGrid grid = grids[gi]
    cdef np.ndarray[np.int32_t, ndim=3] child_indices = grid.child_indices
    cdef np.ndarray[np.int32_t, ndim=1] dimensions = grid.dimensions
    cdef np.ndarray[np.float64_t, ndim=4] fields = grid.fields
    cdef np.ndarray[np.float64_t, ndim=1] leftedges = grid.left_edges
    cdef np.float64_t dx = grid.dx[0]
    cdef np.float64_t child_dx
    cdef np.ndarray[np.float64_t, ndim=1] child_leftedges
    cdef np.float64_t cx, cy, cz
    cdef np.ndarray[np.int32_t, ndim=1] oct
    cdef hilbert_state hs_child
    for oct,hs_child in hs:
        i,j,k = oct[0],oct[1],oct[2]
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
            RecurseOctreeDepthFirst(child_i, child_j, child_k, 2, 2, 2,
                                curpos, ci - grid.offset, hs_child, output, refined, grids)


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
    cdef np.ndarray[np.int32_t, ndim=3] child_indices = grid.child_indices
    cdef np.ndarray[np.int32_t, ndim=1] dimensions = grid.dimensions
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
    cdef np.ndarray[np.int32_t, ndim=1] dimensions = grid.dimensions
    cdef np.ndarray[np.float64_t, ndim=4] fields = grid.fields
    cdef np.ndarray[np.float64_t, ndim=1] leftedges = grid.left_edges
    cdef np.float64_t dx = grid.dx[0]
    cdef np.float64_t child_dx
    cdef np.ndarray[np.float64_t, ndim=1] child_leftedges
    cdef np.float64_t cx, cy, cz
    cdef int cp
    for i_off in range(i_f):
        i = i_off + i_i
        cx = (leftedges[0] + i*dx)
        if i_f > 2: print k, cz
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

