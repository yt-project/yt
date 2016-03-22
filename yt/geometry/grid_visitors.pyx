"""
Grid visitor functions




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free
from yt.utilities.lib.fp_utils cimport iclip
from yt.geometry.grid_container cimport GridTree

cdef class GridVisitor:
    def __cinit__(self):
        self.index = 0
        self.global_index = 0
        self.n_tuples = 0
        self.child_tuples = NULL
        self.ref_factor = 2 #### FIX THIS

    cdef void free_tuples(self) nogil:
        # This wipes out the tuples, which is necessary since they are
        # heap-allocated
        cdef int i
        if self.child_tuples == NULL: return
        for i in range(self.n_tuples):
            free(self.child_tuples[i])
        free(self.child_tuples)
        self.child_tuples = NULL
        self.n_tuples = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void setup_tuples(self, GridTreeNode *grid) nogil:
        # This sets up child-mask tuples.  Rather than a single mask that covers
        # everything, we instead allocate pairs of integers that are start/stop
        # positions for child masks.  This may not be considerably more efficient
        # memory-wise, but it is easier to keep and save when going through
        # multiple grids and selectors.
        cdef int i, j
        cdef np.int64_t si, ei
        cdef GridTreeNode *c
        self.free_tuples()
        self.child_tuples = <int**> malloc(sizeof(int*) * grid.num_children)
        for i in range(grid.num_children):
            c = grid.children[i]
            self.child_tuples[i] = <int *>malloc(sizeof(int) * 6)
            # Now we fill them in
            for j in range(3):
                si = (c.start_index[j] / self.ref_factor) - grid.start_index[j]
                ei = si + c.dims[j]/self.ref_factor - 1
                self.child_tuples[i][j*2+0] = iclip(si, 0, grid.dims[j] - 1)
                self.child_tuples[i][j*2+1] = iclip(ei, 0, grid.dims[j] - 1)
        self.n_tuples = grid.num_children

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.uint8_t check_child_masked(self) nogil:
        # This simply checks if we're inside any of the tuples.  Probably not the
        # most efficient way, but the GVD* passed in has a position affiliated with
        # it, and we can very easily look for that inside here.
        cdef int i, j, k
        cdef int *tup
        for i in range(self.n_tuples):
            # k is if we're inside a given child tuple.  We check each one
            # individually, and invalidate if we're outside.
            k = 1
            tup = self.child_tuples[i]
            for j in range(3):
                # Check if pos is outside in any of the three dimensions
                if self.pos[j] < tup[j*2+0] or self.pos[j] > tup[j*2+1]:
                    k = 0
                    break
            if k == 1: return 1 # Return 1 for child masked
        return 0 # Only return 0 if it doesn't match any of the children

    cdef void visit(self, GridTreeNode *grid, np.uint8_t selected) nogil:
        with gil:
            raise NotImplementedError

cdef class CountGridCells(GridVisitor):
    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void visit(self, GridTreeNode *grid, np.uint8_t selected) nogil:
        # Simply increment for each one, if we've selected it.
        if selected == 0: return
        self.count += 1

cdef class MaskGridCells(GridVisitor):
    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void visit(self, GridTreeNode *grid, np.uint8_t selected) nogil:
        # Set our bitarray -- we're creating a mask -- if we are selected.
        if selected == 0: return
        self.mask[self.global_index] = 1
        # No need to increment anything.

cdef class ICoordsGrids(GridVisitor):
    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void visit(self, GridTreeNode *grid, np.uint8_t selected) nogil:
        # Nice and easy icoord setter.
        if selected == 0: return
        cdef int i
        for i in range(3):
            self.icoords[self.index,i] = self.pos[i] + grid.start_index[i]
        self.index += 1

cdef class IResGrids(GridVisitor):
    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void visit(self, GridTreeNode *grid, np.uint8_t selected) nogil:
        # Fill with the level value.
        if selected == 0: return
        self.ires[self.index] = grid.level
        self.index += 1

cdef class FWidthGrids(GridVisitor):
    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void visit(self, GridTreeNode *grid, np.uint8_t selected) nogil:
        # Fill with our dds.
        if selected == 0: return
        cdef int i
        for i in range(3):
            self.fwidth[self.index,i] = grid.dds[i]
        self.index += 1

cdef class FCoordsGrids(GridVisitor):
    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void visit(self, GridTreeNode *grid, np.uint8_t selected) nogil:
        # Simple cell-centered position filling.
        if selected == 0: return
        cdef int i
        for i in range(3):
            self.fcoords[self.index,i] = grid.left_edge[i] + \
                (0.5 + self.pos[i])*grid.dds[i]
        self.index += 1
