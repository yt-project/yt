"""
Geometry selection routines.

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
from fp_utils cimport fclip
from cython.parallel import prange, parallel, threadid

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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def convert_mask_to_indices(np.ndarray[np.uint8_t, ndim=3, cast=True] mask,
            int count, int transpose = 0):
    cdef int i, j, k, cpos
    cdef np.ndarray[np.int32_t, ndim=2] indices 
    indices = np.zeros((count, 3), dtype='int32')
    cpos = 0
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            for k in range(mask.shape[2]):
                if mask[i,j,k] == 1:
                    if transpose == 1:
                        indices[cpos, 0] = k
                        indices[cpos, 1] = j
                        indices[cpos, 2] = i
                    else:
                        indices[cpos, 0] = i
                        indices[cpos, 1] = j
                        indices[cpos, 2] = k
                    cpos += 1
    return indices

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
    return gridi.astype("bool")

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
    return gridi.astype("bool")

def slice_grids(dobj, np.ndarray[np.float64_t, ndim=2] left_edges,
                      np.ndarray[np.float64_t, ndim=2] right_edges):
    cdef int i, ax
    cdef int ng = left_edges.shape[0]
    cdef np.ndarray[np.int32_t, ndim=1] gridi = np.zeros(ng, dtype='int32')
    ax = dobj.axis
    cdef np.float64_t coord = dobj.coord
    for i in range(ng):
        if right_edges[i, ax] > coord and left_edges[i, ax] <= coord:
            gridi[i] = 1
    return gridi.astype("bool")

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
    return gridi.astype("bool")

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
    return mask.astype("bool")

# Inclined Box
# Sphere

cdef class SelectorObject:

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges):
        cdef int i, n
        cdef int ng = left_edges.shape[0]
        cdef np.ndarray[np.uint8_t, ndim=1] gridi = np.zeros(ng, dtype='uint8')
        cdef np.float64_t LE[3], RE[3]
        with nogil, parallel():
            for n in prange(ng):
                # Call our selector function
                # Check if the sphere is inside the grid
                for i in range(3):
                    LE[i] = left_edges[n, i]
                    RE[i] = right_edges[n, i]
                gridi[n] = self.select_grid(LE, RE)
        return gridi.astype("bool")

    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        return 0
    
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_cells(self, gobj):
        cdef np.ndarray[np.float64_t, ndim=1] odds = gobj.dds
        cdef np.ndarray[np.float64_t, ndim=1] left_edge = gobj.LeftEdge
        cdef np.ndarray[np.float64_t, ndim=1] right_edge = gobj.RightEdge
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef np.float64_t dds[3], pos[3]
        cdef int i, j, k, nv[3]
        child_mask = gobj.child_mask
        for i in range(3):
            nv[i] = gobj.ActiveDimensions[i]
            dds[i] = odds[i]
        cdef int count = 0
        with nogil:
            pos[0] = left_edge[0] + dds[0] * 0.5
            for i in range(nv[0]):
                pos[1] = left_edge[1] + dds[1] * 0.5
                for j in range(nv[1]):
                    pos[2] = left_edge[2] + dds[2] * 0.5
                    for k in range(nv[2]):
                        if child_mask[i,j,k] == 1:
                            count += self.select_cell(pos, dds)
                        pos[2] += dds[1]
                    pos[1] += dds[1]
                pos[0] += dds[0]
        return count

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        child_mask = gobj.child_mask
        cdef np.ndarray[np.uint8_t, ndim=3] mask 
        cdef int nv[3]
        cdef np.ndarray[np.float64_t, ndim=1] odds = gobj.dds
        cdef np.ndarray[np.float64_t, ndim=1] left_edge = gobj.LeftEdge
        cdef np.ndarray[np.float64_t, ndim=1] right_edge = gobj.RightEdge
        cdef int i, j, k
        cdef np.float64_t dds[3], pos[3]
        for i in range(3):
            dds[i] = odds[i]
            nv[i] = gobj.ActiveDimensions[i]
        mask = np.zeros(gobj.ActiveDimensions, dtype='uint8')
        cdef int temp
        with nogil:
            pos[0] = left_edge[0] + dds[0] * 0.5
            for i in range(nv[0]):
                pos[1] = left_edge[1] + dds[1] * 0.5
                for j in range(nv[1]):
                    pos[2] = left_edge[2] + dds[2] * 0.5
                    for k in range(nv[2]):
                        if child_mask[i,j,k] == 1:
                            mask[i,j,k] = self.select_cell(pos, dds)
                        pos[2] += dds[1]
                    pos[1] += dds[1]
                pos[0] += dds[0]
        return mask.astype("bool")

cdef class SphereSelector(SelectorObject):
    cdef np.float64_t radius2
    cdef np.float64_t center[3]

    def __init__(self, dobj):
        for i in range(3):
            self.center[i] = dobj.center[i]
        self.radius2 = dobj.radius * dobj.radius

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef np.float64_t box_center, relcenter, closest, dist, edge
        cdef int id
        if (left_edge[0] <= self.center[0] <= right_edge[0] and
            left_edge[1] <= self.center[1] <= right_edge[1] and
            left_edge[2] <= self.center[2] <= right_edge[2]):
            return 1
        # http://www.gamedev.net/topic/335465-is-this-the-simplest-sphere-aabb-collision-test/
        dist = 0
        for i in range(3):
            box_center = (right_edge[i] + left_edge[i])/2.0
            relcenter = self.center[i] - box_center
            edge = right_edge[i] - left_edge[i]
            closest = relcenter - fclip(relcenter, -edge/2.0, edge/2.0)
            dist += closest * closest
        if dist < self.radius2: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef np.float64_t dist2, temp
        cdef int i
        dist2 = 0
        for i in range(3):
            temp = (pos[i] - self.center[i])
            dist2 += temp * temp
        if dist2 <= self.radius2: return 1
        return 0

sphere_selector = SphereSelector

cdef class RegionSelector(SelectorObject):
    cdef np.float64_t left_edge[3]
    cdef np.float64_t right_edge[3]
    cdef np.float64_t dx_pad

    def __init__(self, dobj):
        cdef int i
        self.dx_pad =dobj._dx_pad
        for i in range(3):
            self.left_edge[i] = dobj.left_edge[i]
            self.right_edge[i] = dobj.right_edge[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        for i in range(3):
            if left_edge[i] >= self.right_edge[i]: return 0
            if right_edge[i] <= self.left_edge[i]: return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef np.float64_t dxp
        for i in range(3):
            dxp = self.dx_pad * dds[i]
            if pos[i] - dxp >= self.right_edge[i]: return 0
            if pos[i] + dxp <= self.left_edge[i]: return 0
        return 1

region_selector = RegionSelector

cdef class DiskSelector(SelectorObject):
    cdef np.float64_t norm_vec[3]
    cdef np.float64_t center[3]
    cdef np.float64_t d
    cdef np.float64_t radius2
    cdef np.float64_t height

    def __init__(self, dobj):
        cdef int i
        for i in range(3):
            self.norm_vec[i] = dobj._norm_vec[i]
            self.center[i] = dobj.center[i]
        self.d = dobj._d
        self.radius2 = dobj._radius * dobj._radius
        self.height = dobj._height

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef np.float64_t *arr[2]
        cdef np.float64_t pos[3], H, D, R2, temp
        cdef int i, j, k, n
        arr[0] = left_edge
        arr[1] = right_edge
        cdef int cond[2]
        cond[0] = cond[1] = 0
        for i in range(2):
            pos[0] = arr[i][0]
            for j in range(2):
                pos[1] = arr[j][1]
                for k in range(2):
                    pos[2] = arr[k][2]
                    H = D = 0
                    for n in range(3):
                        H += (pos[n] * self.norm_vec[n])
                        temp = (pos[n] - self.center[n])
                        D += temp*temp
                    H += self.d
                    R2 = (D - H*H)
                    if fabs(H) < self.height: cond[0] = 1
                    if R2 < self.radius2: cond[1] = 1
        # A moment of explanation:
        #    We want our height to be less than the height AND our radius2 to be
        #    less than radius2, so we set cond[0] equal to 1 if any corners
        #    match that criteria.
        # Note that we OVERSELECT grids, as we are selecting anything within
        # the height and within the radius, which is kind of a funny thing.
        # Cell selection takes care of the rest.
        if cond[0] == 1 and cond[1] == 1:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef np.float64_t h, d, r2, temp
        cdef int i
        h = d = 0
        for i in range(3):
            h += pos[i] * self.norm_vec[i]
            temp = pos[i] - self.center[i]
            d += temp*temp
        h += self.d
        r2 = (d - h*h)
        if fabs(h) <= self.height and r2 <= self.radius2: return 1
        return 0

disk_selector = DiskSelector

# Ellipse
