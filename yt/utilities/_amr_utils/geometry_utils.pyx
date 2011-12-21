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
from stdlib cimport malloc, free

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

# These routines are separated into a couple different categories:
#
#   * Routines for identifying intersections of an object with a bounding box
#   * Routines for identifying cells/points inside a bounding box that
#     intersect with an object
#   * Routines that speed up some type of geometric calculation

# First, bounding box / object intersection routines.
# These all respect the interface "dobj" and a set of left_edges, right_edges,
# sometimes also accepting level and mask information.

def ortho_ray_grids(dobj, np.ndarray[np.float64_t, ndim=2] left_edges,
                          np.ndarray[np.float64_t, ndim=2] right_edges):
    cdef int i
    cdef int ng = left_edges.shape[0]
    cdef int px_ax = dobj.px_ax
    cdef int py_ax = dobj.py_ax
    cdef np.float64_t px = dobj.px
    cdef np.float64_t py = dobj.py
    cdef np.ndarray[np.int32_t, ndim=1] gridi = np.zeros(ng, dtype='int32_t')
    for i in range(ng):
        if (    (px >= left_edges[i, px])
            and (px < right_edges[i, px])
            and (py >= left_edges[i, py])
            and (py < right_edges[i, py])):
            gridi[i] = 1
    return gridi

def ray_grids(dobj, np.ndarray[np.float64_t, ndim=2] left_edges,
                    np.ndarray[np.float64_t, ndim=2] right_edges):
    cdef int i, ax
    cdef int i1, i2
    cdef int ng = left_edges.shape[0]
    cdef np.ndarray[np.int32_t, ndim=1] gridi = np.zeros(ng, dtype='int32')
    cdef np.float64_t vs[3], t, p0[3], p1[3], v[3]
    for i in range(3):
        p0[i] = dobj.start_point[i]
        p1[i] = dobj.end_point[i]
        v[i] = dobj.vec[i]
    # We check first to see if at any point, the ray intersects a grid face
    for gi in range(ng):
        for ax in range(3):
            i1 = (ax+1) % 3
            i2 = (ax+2) % 3
            t = (left_edges[gi,ax] - p0[ax])/v[ax]
            for i in range(3):
                vs[i] = t * v[i] + p0[i]
            if left_edges[gi,i1] <= vs[i1] and \
               right_edges[gi,i1] >= vs[i1] and \
               left_edges[gi,i2] <= vs[i2] and \
               right_edges[gi,i2] >= vs[i2]:
                gridi[gi] = 1
                break
            t = (right_edges[gi,ax] - p0[ax])/v[ax]
            for i in range(3):
                vs[i] = t * v[i] + p0[i]
            if left_edges[gi,i1] <= vs[i1] and \
               right_edges[gi,i1] >= vs[i1] and \
               left_edges[gi,i2] <= vs[i2] and \
               right_edges[gi,i2] >= vs[i2]:
                gridi[gi] = 1
                break
        if gridi[gi] == 1: continue
        # if the point is fully enclosed, we count the grid
        if left_edges[gi,0] <= p0[0] and \
           right_edges[gi,0] >= p0[0] and \
           left_edges[gi,1] <= p0[1] and \
           right_edges[gi,1] >= p0[1] and \
           left_edges[gi,2] <= p0[2] and \
           right_edges[gi,2] >= p0[2]:
            gridi[gi] = 1
            continue
        if left_edges[gi,0] <= p1[0] and \
           right_edges[gi,0] >= p1[0] and \
           left_edges[gi,1] <= p1[1] and \
           right_edges[gi,1] >= p1[1] and \
           left_edges[gi,2] <= p1[2] and \
           right_edges[gi,2] >= p1[2]:
            gridi[gi] = 1
            continue
    return gridi

def slice_grids(dobj, np.ndarray[np.float64_t, ndim=2] left_edges,
                      np.ndarray[np.float64_t, ndim=2] right_edges):
    cdef int i, ax
    cdef int ng = left_edges.shape[0]
    cdef np.ndarray[np.int32_t, ndim=1] gridi = np.zeros(ng, dtype='int32')
    ax = dobj.axis
    cdef np.float64_t coord = dobj.coord
    for i in range(ng):
        if left_edges[i, ax] <= coord and \
           right_edges[i, ax] > coord:
            gridi[i] = 1
    return gridi

def cutting_plane_grids(dobj, np.ndarray[np.float64_t, ndim=2] left_edges,
                        np.ndarray[np.float64_t, ndim=2] right_edges):
    cdef int i
    cdef int ng = left_edges.shape[0]
    cdef np.ndarray[np.int32_t, ndim=1] gridi = np.zeros(ng, dtype='int32')
    cdef np.float64_t *arr[2]
    arr[0] = <np.float64_t *> left_edges.data
    arr[1] = <np.float64_t *> right_edges.data
    cdef np.float64_t x, y, z
    cdef np.float64_t norm_vec[3]
    cdef np.float64_t d = dobj._d # offset to center
    cdef np.float64_t gd # offset to center
    cdef np.int64_t all_under, all_over
    for i in range(3):
        norm_vec[i] = dobj._norm_vec[i]
    for i in range(ng):
        all_under = 1
        all_over = 1
        # Check each corner
        for xi in range(2):
            x = arr[xi][i * 3 + 0]
            for yi in range(2):
                y = arr[yi][i * 3 + 1]
                for zi in range(2):
                    z = arr[zi][i * 3 + 2]
                    gd = ( x*norm_vec[0]
                         + y*norm_vec[1]
                         + z*norm_vec[2]) + d
                    if gd <= 0: all_over = 0
                    if gd >= 0: all_under = 0
        if not (all_over == 1 or all_under == 1):
            gridi[i] = 1
    return gridi

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline int cutting_plane_cell(
                        np.float64_t x, np.float64_t y, np.float64_t z,
                        np.float64_t norm_vec[3], np.float64_t d,
                        np.float64_t dist):
    cdef np.float64_t cd = x*norm_vec[0] + y*norm_vec[1] + z*norm_vec[2] + d
    if fabs(cd) <= dist: return 1
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cutting_plane_cells(dobj, gobj):
    cdef np.ndarray[np.int32_t, ndim=3] mask 
    cdef np.ndarray[np.float64_t, ndim=1] left_edge = gobj.LeftEdge
    cdef np.ndarray[np.float64_t, ndim=1] dds = gobj.dds
    cdef int i, j, k
    cdef np.float64_t x, y, z, dist
    cdef np.float64_t norm_vec[3]
    cdef np.float64_t d = dobj._d

    mask = np.zeros(gobj.ActiveDimensions, dtype='int32')
    for i in range(3): norm_vec[i] = dobj._norm_vec[i]
    dist = 0.5*(dds[0]*dds[0] + dds[1]*dds[1] + dds[2]*dds[2])**0.5
    x = left_edge[0] + dds[0] * 0.5
    for i in range(mask.shape[0]):
        y = left_edge[1] + dds[1] * 0.5
        for j in range(mask.shape[1]):
            z = left_edge[2] + dds[2] * 0.5
            for k in range(mask.shape[2]):
                mask[i,j,k] = cutting_plane_cell(x, y, z, norm_vec, d, dist)
                z += dds[1]
            y += dds[1]
        x += dds[0]
    return mask
                
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
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

