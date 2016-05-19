"""
Utilities for unstructured and semi-structured meshes



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
from libc.stdlib cimport malloc, free, abs
from yt.utilities.lib.fp_utils cimport imax, fmax, imin, fmin, iclip, fclip, i64clip

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def fill_fcoords(np.ndarray[np.float64_t, ndim=2] coords,
                 np.ndarray[np.int64_t, ndim=2] indices,
                 int offset = 0):
    cdef np.ndarray[np.float64_t, ndim=2] fcoords
    cdef int nc = indices.shape[0]
    cdef int nv = indices.shape[1]
    cdef np.float64_t pos[3]
    cdef int i, j, k
    fcoords = np.empty((nc, 3), dtype="float64")
    for i in range(nc):
        for j in range(3):
            pos[j] = 0.0
        for j in range(nv):
            for k in range(3):
                pos[k] += coords[indices[i, j] - offset, k]
        for j in range(3):
            pos[j] /= nv
            fcoords[i, j] = pos[j]
    return fcoords

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def fill_fwidths(np.ndarray[np.float64_t, ndim=2] coords,
                 np.ndarray[np.int64_t, ndim=2] indices,
                 int offset = 0):
    cdef np.ndarray[np.float64_t, ndim=2] fwidths
    cdef int nc = indices.shape[0]
    cdef int nv = indices.shape[1]
    if nv != 8:
        raise NotImplementedError
    cdef np.float64_t LE[3]
    cdef np.float64_t RE[3]
    cdef int i, j, k
    cdef np.float64_t pos
    fwidths = np.empty((nc, 3), dtype="float64")
    for i in range(nc):
        for j in range(3):
            LE[j] = 1e60
            RE[j] = -1e60
        for j in range(nv):
            for k in range(3):
                pos = coords[indices[i, j] - offset, k]
                LE[k] = fmin(pos, LE[k])
                RE[k] = fmax(pos, RE[k])
        for j in range(3):
            fwidths[i, j] = RE[j] - LE[j]
    return fwidths

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def smallest_fwidth(np.ndarray[np.float64_t, ndim=2] coords,
                    np.ndarray[np.int64_t, ndim=2] indices,
                    int offset = 0):
    cdef np.float64_t fwidth = 1e60, pos
    cdef int nc = indices.shape[0]
    cdef int nv = indices.shape[1]
    if nv != 8:
        raise NotImplementedError
    cdef np.float64_t LE[3]
    cdef np.float64_t RE[3]
    cdef int i, j, k
    for i in range(nc):
        for j in range(3):
            LE[j] = 1e60
            RE[j] = -1e60
        for j in range(nv):
            for k in range(3):
                pos = coords[indices[i, j] - offset, k]
                LE[k] = fmin(pos, LE[k])
                RE[k] = fmax(pos, RE[k])
        for j in range(3):
            fwidth = fmin(fwidth, RE[j] - LE[j])
    return fwidth

