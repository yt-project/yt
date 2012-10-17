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
from libc.stdlib cimport malloc, free
from fp_utils cimport fclip

cdef extern from "math.h":
    double exp(double x) nogil
    float expf(float x) nogil
    long double expl(long double x) nogil
    double floor(double x) nogil
    double ceil(double x) nogil
    double fmod(double x, double y) nogil
    double log2(double x) nogil
    long int lrint(double x) nogil
    double fabs(double x) nogil

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def get_box_grids_level(np.ndarray[np.float64_t, ndim=1] left_edge,
                        np.ndarray[np.float64_t, ndim=1] right_edge,
                        int level,
                        np.ndarray[np.float64_t, ndim=2] left_edges,
                        np.ndarray[np.float64_t, ndim=2] right_edges,
                        np.ndarray[np.int32_t, ndim=2] levels,
                        np.ndarray[np.int32_t, ndim=1] mask,
                        int min_index = 0):
    cdef int i, n
    cdef int nx = left_edges.shape[0]
    cdef int inside 
    for i in range(nx):
        if i < min_index or levels[i,0] != level:
            mask[i] = 0
            continue
        inside = 1
        for n in range(3):
            if left_edge[n] >= right_edges[i,n] or \
               right_edge[n] <= left_edges[i,n]:
                inside = 0
                break
        if inside == 1: mask[i] = 1
        else: mask[i] = 0

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def get_box_grids_below_level(
                        np.ndarray[np.float64_t, ndim=1] left_edge,
                        np.ndarray[np.float64_t, ndim=1] right_edge,
                        int level,
                        np.ndarray[np.float64_t, ndim=2] left_edges,
                        np.ndarray[np.float64_t, ndim=2] right_edges,
                        np.ndarray[np.int32_t, ndim=2] levels,
                        np.ndarray[np.int32_t, ndim=1] mask):
    cdef int i, n
    cdef int nx = left_edges.shape[0]
    cdef int inside 
    for i in range(nx):
        mask[i] = 0
        if levels[i,0] <= level:
            inside = 1
            for n in range(3):
                if left_edge[n] >= right_edges[i,n] or \
                   right_edge[n] <= left_edges[i,n]:
                    inside = 0
                    break
            if inside == 1: mask[i] = 1

# Finally, miscellaneous routines.

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def find_values_at_point(np.ndarray[np.float64_t, ndim=1] point,
                         np.ndarray[np.float64_t, ndim=2] left_edges,
                         np.ndarray[np.float64_t, ndim=2] right_edges,
                         np.ndarray[np.int32_t, ndim=2] dimensions,
                         field_names, grid_objects):
    # This iterates in order, first to last, and then returns with the first
    # one in which the point is located; this means if you order from highest
    # level to lowest, you will find the correct grid without consulting child
    # masking.  Note also that we will do a few relatively slow operations on
    # strings and whatnot, but they should not be terribly slow.
    cdef int ind[3], gi, fi
    cdef int nf = len(field_names)
    cdef np.float64_t dds
    cdef np.ndarray[np.float64_t, ndim=3] field
    cdef np.ndarray[np.float64_t, ndim=1] rv = np.zeros(nf, dtype='float64')
    for gi in range(left_edges.shape[0]):
        if not ((left_edges[gi,0] < point[0] < right_edges[gi,0])
            and (left_edges[gi,1] < point[1] < right_edges[gi,1])
            and (left_edges[gi,2] < point[2] < right_edges[gi,2])):
            continue
        # We found our grid!
        for fi in range(3):
            dds = ((right_edges[gi,fi] - left_edges[gi,fi])/
                   (<np.float64_t> dimensions[gi,fi]))
            ind[fi] = <int> ((point[fi] - left_edges[gi,fi])/dds)
        grid = grid_objects[gi]
        for fi in range(nf):
            field = grid[field_names[fi]]
            rv[fi] = field[ind[0], ind[1], ind[2]]
        return rv
    raise KeyError

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def obtain_rvec(data):
    # This is just to let the pointers exist and whatnot.  We can't cdef them
    # inside conditionals.
    cdef np.ndarray[np.float64_t, ndim=1] xf
    cdef np.ndarray[np.float64_t, ndim=1] yf
    cdef np.ndarray[np.float64_t, ndim=1] zf
    cdef np.ndarray[np.float64_t, ndim=2] rf
    cdef np.ndarray[np.float64_t, ndim=3] xg
    cdef np.ndarray[np.float64_t, ndim=3] yg
    cdef np.ndarray[np.float64_t, ndim=3] zg
    cdef np.ndarray[np.float64_t, ndim=4] rg
    cdef np.float64_t c[3]
    cdef int i, j, k
    center = data.get_field_parameter("center")
    c[0] = center[0]; c[1] = center[1]; c[2] = center[2]
    if len(data['x'].shape) == 1:
        # One dimensional data
        xf = data['x']
        yf = data['y']
        zf = data['z']
        rf = np.empty((3, xf.shape[0]), 'float64')
        for i in range(xf.shape[0]):
            rf[0, i] = xf[i] - c[0]
            rf[1, i] = yf[i] - c[1]
            rf[2, i] = zf[i] - c[2]
        return rf
    else:
        # Three dimensional data
        xg = data['x']
        yg = data['y']
        zg = data['z']
        rg = np.empty((3, xg.shape[0], xg.shape[1], xg.shape[2]), 'float64')
        for i in range(xg.shape[0]):
            for j in range(xg.shape[1]):
                for k in range(xg.shape[2]):
                    rg[0,i,j,k] = xg[i,j,k] - c[0]
                    rg[1,i,j,k] = yg[i,j,k] - c[1]
                    rg[2,i,j,k] = zg[i,j,k] - c[2]
        return rg

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t graycode(np.int64_t x):
    return x^(x>>1)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t igraycode(np.int64_t x):
    cdef np.int64_t i, j
    if x == 0:
        return x
    m = <np.int64_t> ceil(log2(x)) + 1
    i, j = x, 1
    while j < m:
        i = i ^ (x>>j)
        j += 1
    return i

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t direction(np.int64_t x, np.int64_t n):
    #assert x < 2**n
    if x == 0:
        return 0
    elif x%2 == 0:
        return tsb(x-1, n)%n
    else:
        return tsb(x, n)%n

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t tsb(np.int64_t x, np.int64_t width):
    #assert x < 2**width
    cdef np.int64_t i = 0
    while x&1 and i <= width:
        x = x >> 1
        i += 1
    return i

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t bitrange(np.int64_t x, np.int64_t width,
                         np.int64_t start, np.int64_t end):
    return x >> (width-end) & ((2**(end-start))-1)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t rrot(np.int64_t x, np.int64_t i, np.int64_t width):
    i = i%width
    x = (x>>i) | (x<<width-i)
    return x&(2**width-1)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t lrot(np.int64_t x, np.int64_t i, np.int64_t width):
    i = i%width
    x = (x<<i) | (x>>width-i)
    return x&(2**width-1)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t transform(np.int64_t entry, np.int64_t direction,
                          np.int64_t width, np.int64_t x):
    return rrot((x^entry), direction + 1, width)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t entry(np.int64_t x):
    if x == 0: return 0
    return graycode(2*((x-1)/2))

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t setbit(np.int64_t x, np.int64_t w, np.int64_t i, np.int64_t b):
    if b == 1:
        return x | 2**(w-i-1)
    elif b == 0:
        return x & ~2**(w-i-1)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t point_to_hilbert(int order, np.int64_t p[3]):
    cdef np.int64_t h, e, d, l, b, w, i, x
    h = e = d = 0
    for i in range(order):
        l = 0
        for x in range(3):
            b = bitrange(p[3-x-1], order, i, i+1)
            l |= (b<<x)
        l = transform(e, d, 3, l)
        w = igraycode(l)
        e = e ^ lrot(entry(w), d+1, 3)
        d = (d + direction(w, 3) + 1)%3
        h = (h<<3)|w
    return h

#def hilbert_point(dimension, order, h):
#    """
#        Convert an index on the Hilbert curve of the specified dimension and
#        order to a set of point coordinates.
#    """
#    #    The bit widths in this function are:
#    #        p[*]  - order
#    #        h     - order*dimension
#    #        l     - dimension
#    #        e     - dimension
#    hwidth = order*dimension
#    e, d = 0, 0
#    p = [0]*dimension
#    for i in range(order):
#        w = utils.bitrange(h, hwidth, i*dimension, i*dimension+dimension)
#        l = utils.graycode(w)
#        l = itransform(e, d, dimension, l)
#        for j in range(dimension):
#            b = utils.bitrange(l, dimension, j, j+1)
#            p[j] = utils.setbit(p[j], order, i, b)
#        e = e ^ utils.lrot(entry(w), d+1, dimension)
#        d = (d + direction(w, dimension) + 1)%dimension
#    return p

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void hilbert_to_point(int order, np.int64_t h, np.int64_t *p):
    cdef np.int64_t hwidth, e, d, w, l, b
    cdef int i, j
    hwidth = 3 * order
    e = d = p[0] = p[1] = p[2] = 0
    for i in range(order):
        w = bitrange(h, hwidth, i*3, i*3+3)
        l = graycode(w)
        l = lrot(l, d +1, 3)^e
        for j in range(3):
            b = bitrange(l, 3, j, j+1)
            p[j] = setbit(p[j], order, i, b)
        e = e ^ lrot(entry(w), d+1, 3)
        d = (d + direction(w, 3) + 1)%3

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def get_hilbert_indices(int order, np.ndarray[np.int64_t, ndim=2] left_index):
    # This is inspired by the scurve package by user cortesi on GH.
    cdef int i
    cdef np.int64_t p[3]
    cdef np.ndarray[np.int64_t, ndim=1] hilbert_indices
    hilbert_indices = np.zeros(left_index.shape[0], 'int64')
    for i in range(left_index.shape[0]):
        p[0] = left_index[i, 0]
        p[1] = left_index[i, 1]
        p[2] = left_index[i, 2]
        hilbert_indices[i] = point_to_hilbert(order, p)
    return hilbert_indices

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def get_hilbert_points(int order, np.ndarray[np.int64_t, ndim=1] indices):
    # This is inspired by the scurve package by user cortesi on GH.
    cdef int i, j
    cdef np.int64_t p[3]
    cdef np.ndarray[np.int64_t, ndim=2] positions
    positions = np.zeros((indices.shape[0], 3), 'int64')
    for i in range(indices.shape[0]):
        hilbert_to_point(order, indices[i], p)
        for j in range(3):
            positions[i, j] = p[j]
    return positions


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def obtain_rv_vec(data):
    # This is just to let the pointers exist and whatnot.  We can't cdef them
    # inside conditionals.
    cdef np.ndarray[np.float64_t, ndim=1] vxf
    cdef np.ndarray[np.float64_t, ndim=1] vyf
    cdef np.ndarray[np.float64_t, ndim=1] vzf
    cdef np.ndarray[np.float64_t, ndim=2] rvf
    cdef np.ndarray[np.float64_t, ndim=3] vxg
    cdef np.ndarray[np.float64_t, ndim=3] vyg
    cdef np.ndarray[np.float64_t, ndim=3] vzg
    cdef np.ndarray[np.float64_t, ndim=4] rvg
    cdef np.float64_t bv[3]
    cdef int i, j, k
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = np.zeros(3)
    bv[0] = bulk_velocity[0]; bv[1] = bulk_velocity[1]; bv[2] = bulk_velocity[2]
    if len(data['x-velocity'].shape) == 1:
        # One dimensional data
        vxf = data['x-velocity'].astype("float64")
        vyf = data['y-velocity'].astype("float64")
        vzf = data['z-velocity'].astype("float64")
        rvf = np.empty((3, vxf.shape[0]), 'float64')
        for i in range(vxf.shape[0]):
            rvf[0, i] = vxf[i] - bv[0]
            rvf[1, i] = vyf[i] - bv[1]
            rvf[2, i] = vzf[i] - bv[2]
        return rvf
    else:
        # Three dimensional data
        vxg = data['x-velocity'].astype("float64")
        vyg = data['y-velocity'].astype("float64")
        vzg = data['z-velocity'].astype("float64")
        rvg = np.empty((3, vxg.shape[0], vxg.shape[1], vxg.shape[2]), 'float64')
        for i in range(vxg.shape[0]):
            for j in range(vxg.shape[1]):
                for k in range(vxg.shape[2]):
                    rvg[0,i,j,k] = vxg[i,j,k] - bv[0]
                    rvg[1,i,j,k] = vyg[i,j,k] - bv[1]
                    rvg[2,i,j,k] = vzg[i,j,k] - bv[2]
        return rvg
