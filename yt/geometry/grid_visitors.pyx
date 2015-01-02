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

cdef void free_tuples(GridVisitorData *data) nogil:
    cdef int i
    if data.child_tuples == NULL: return
    for i in range(data.n_tuples):
        free(data.child_tuples[i])
    free(data.child_tuples)
    data.child_tuples = NULL
    data.n_tuples = 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void setup_tuples(GridVisitorData *data) nogil:
    cdef int i, j, k
    cdef np.int64_t si, ei
    cdef GridTreeNode *g, *c
    free_tuples(data)
    g = data.grid
    data.child_tuples = <int**> malloc(sizeof(int*) * g.num_children)
    for i in range(g.num_children):
        c = g.children[i]
        data.child_tuples[i] = <int *>malloc(sizeof(int) * 6)
        # Now we fill them in
        for j in range(3):
            si = (c.start_index[j] / data.ref_factor) - g.start_index[j]
            ei = si + c.dims[j]/data.ref_factor - 1
            data.child_tuples[i][j*2+0] = iclip(si, 0, g.dims[j] - 1)
            data.child_tuples[i][j*2+1] = iclip(ei, 0, g.dims[j] - 1)
    data.n_tuples = g.num_children

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.uint8_t check_child_masked(GridVisitorData *data) nogil:
    cdef int i, j, k
    cdef int *tup
    for i in range(data.n_tuples):
        k = 1
        tup = data.child_tuples[i]
        for j in range(3):
            # Check if pos is outside in any of the three dimensions
            if data.pos[j] < tup[j*2+0] or data.pos[j] > tup[j*2+1]:
                k = 0
                break
        if k == 1: return 1 # Return 1 for child masked
    return 0 # Only return 0 if it doesn't match any of the children

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void count_cells(GridVisitorData *data, np.uint8_t selected,
                      np.uint8_t child_masked) nogil:
    if selected == 0 or child_masked == 1: return
    cdef np.uint64_t *count = <np.uint64_t*> data.array
    count[0] += 1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void icoords_cells(GridVisitorData *data, np.uint8_t selected,
                        np.uint8_t child_masked) nogil:
    cdef int i
    if selected == 0 or child_masked == 1: return
    cdef np.int64_t *icoords = <np.int64_t*> data.array 
    for i in range(3):
        icoords[data.index * 3 + i] = data.pos[i] + data.grid.start_index[i]
    data.index += 1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void ires_cells(GridVisitorData *data, np.uint8_t selected,
                     np.uint8_t child_masked) nogil:
    cdef int i
    if selected == 0 or child_masked == 1: return
    cdef np.int64_t *ires = <np.int64_t*> data.array
    ires[data.index] = data.grid.level
    data.index += 1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void fwidth_cells(GridVisitorData *data, np.uint8_t selected,
                        np.uint8_t child_masked) nogil:
    cdef int i
    if selected == 0 or child_masked == 1: return
    cdef np.float64_t *fwidth = <np.float64_t*> data.array 
    for i in range(3):
        fwidth[data.index * 3 + i] = data.grid.dds[i]
    data.index += 1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void fcoords_cells(GridVisitorData *data, np.uint8_t selected,
                        np.uint8_t child_masked) nogil:
    cdef int i
    if selected == 0 or child_masked == 1: return
    cdef np.float64_t *fcoords = <np.float64_t*> data.array 
    for i in range(3):
        fcoords[data.index * 3 + i] = data.grid.left_edge[i] + \
            (0.5 + data.pos[i])*data.grid.dds[i]
    data.index += 1
