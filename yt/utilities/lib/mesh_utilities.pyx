
# distutils: libraries = STD_LIBS
"""
Utilities for unstructured and semi-structured meshes



"""


import numpy as np

cimport cython
cimport numpy as np
from libc.stdlib cimport abs, free, malloc

from yt.utilities.lib.fp_utils cimport fclip, fmax, fmin, i64clip, iclip, imax, imin

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
