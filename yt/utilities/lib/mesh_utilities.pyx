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
from yt.units.yt_array import YTArray

cdef extern from "platform_dep.h":
    double rint(double x)

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

def fill_fcoords_extruded(np.ndarray[np.float64_t, ndim=3] grid_xy,
                          np.ndarray[np.float64_t, ndim=1] grid_z):
    cdef np.ndarray[np.float64_t, ndim=2] fcoords
    # This is for barycenters
    nv = 8
    cdef int nx, ny, nz, ind
    cdef np.float64_t pos[3]
    nx = grid_xy.shape[0]
    ny = grid_xy.shape[1]
    nz = grid_z.shape[0]
    nc = (nx - 1) * (ny - 1) * (nz - 1)
    fcoords = np.empty((nc, 3), dtype="float64")
    ind = 0
    for i in range(nx - 1):
        for j in range(ny - 1):
            pos[0] = (grid_xy[i,j,0] + grid_xy[i+1,j,0] +
                      grid_xy[i+1,j,0] + grid_xy[i+1,j+1,0])/4.0
            pos[1] = (grid_xy[i,j,1] + grid_xy[i+1,j,1] +
                      grid_xy[i+1,j,1] + grid_xy[i+1,j+1,1])/4.0
            for k in range(nz - 1):
                pos[2] = (grid_z[k] + grid_z[k+1])/2.0
                fcoords[ind, 0] = pos[0]
                fcoords[ind, 1] = pos[1]
                fcoords[ind, 2] = pos[2]
                ind += 1
    return fcoords

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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def clamp_edges(np.float64_t[:] edge, np.float64_t[:] pleft,
                np.float64_t[:] pdx):
    """Clamp edge to pleft + pdx*n where n is the closest integer

    Note that edge is modified in-place.
    """
    cdef np.float64_t start_index
    cdef np.float64_t integer_index
    cdef np.intp_t shape = edge.shape[0]
    for i in range(shape):
        start_index = (edge[i] - pleft[i]) / pdx[i]
        integer_index = rint(start_index)
        edge[i] = integer_index * pdx[i] + pleft[i]
