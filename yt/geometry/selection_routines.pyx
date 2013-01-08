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
from fp_utils cimport fclip, iclip
from cython.parallel import prange, parallel, threadid
from selection_routines cimport SelectorObject
from oct_container cimport OctreeContainer, OctAllocationContainer, Oct
#from geometry_utils cimport point_to_hilbert
from yt.utilities.lib.grid_traversal cimport \
    VolumeContainer, sample_function, walk_volume

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

ctypedef fused anyfloat:
    np.float32_t
    np.float64_t

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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef _mask_fill(np.ndarray[np.float64_t, ndim=1] out,
                np.int64_t offset,
                np.ndarray[np.uint8_t, ndim=3, cast=True] mask,
                np.ndarray[anyfloat, ndim=3] vals):
    cdef np.int64_t count = 0
    cdef int i, j, k
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            for k in range(mask.shape[2]):
                if mask[i,j,k] == 1:
                    out[offset + count] = vals[i,j,k]
                    count += 1
    return count

def mask_fill(np.ndarray[np.float64_t, ndim=1] out,
              np.int64_t offset,
              np.ndarray[np.uint8_t, ndim=3, cast=True] mask,
              np.ndarray vals):
    if vals.dtype == np.float32:
        return _mask_fill[np.float32_t](out, offset, mask, vals)
    elif vals.dtype == np.float64:
        return _mask_fill[np.float64_t](out, offset, mask, vals)
    else:
        raise RuntimeError

# Inclined Box

cdef class SelectorObject:

    def __cinit__(self, dobj):
        self.min_level = getattr(dobj, "min_level", 0)
        self.max_level = getattr(dobj, "max_level", 99)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        cdef int i, n
        cdef int ng = left_edges.shape[0]
        cdef np.ndarray[np.uint8_t, ndim=1] gridi = np.zeros(ng, dtype='uint8')
        cdef np.float64_t LE[3], RE[3]
        with nogil:
            for n in range(ng):
                # Call our selector function
                # Check if the sphere is inside the grid
                for i in range(3):
                    LE[i] = left_edges[n, i]
                    RE[i] = right_edges[n, i]
                gridi[n] = self.select_grid(LE, RE, levels[n,0])
        return gridi.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_octs(self, OctreeContainer octree):
        cdef int i, j, k, n
        cdef np.ndarray[np.uint8_t, ndim=2] mask = np.zeros((octree.nocts, 8), dtype='uint8')
        cdef np.float64_t pos[3], dds[3]
        # This dds is the oct-width
        for i in range(3):
            dds[i] = (octree.DRE[i] - octree.DLE[i]) / octree.nn[i]
        # Pos is the center of the octs
        pos[0] = octree.DLE[0] + dds[0]/2.0
        for i in range(octree.nn[0]):
            pos[1] = octree.DLE[1] + dds[1]/2.0
            for j in range(octree.nn[1]):
                pos[2] = octree.DLE[2] + dds[2]/2.0
                for k in range(octree.nn[2]):
                    if octree.root_mesh[i][j][k] == NULL: continue
                    self.recursively_select_octs(
                        octree.root_mesh[i][j][k],
                        pos, dds, mask, 0) # last is level
                    pos[2] += dds[2]
                pos[1] += dds[2]
            pos[0] += dds[2]
        return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void recursively_select_octs(self, Oct *root,
                        np.float64_t pos[3], np.float64_t dds[3],
                        np.ndarray[np.uint8_t, ndim=2] mask,
                        int level = 0):
        cdef np.float64_t LE[3], RE[3], sdds[3], spos[3]
        cdef int i, j, k, res, ii
        cdef Oct *ch
        # Remember that pos is the *center* of the oct, and dds is the oct
        # width.  So to get to the edges, we add/subtract half of dds.
        for i in range(3):
            # sdds is the cell width
            sdds[i] = dds[i]/2.0
            LE[i] = pos[i] - dds[i]/2.0
            RE[i] = pos[i] + dds[i]/2.0
        #print LE[0], RE[0], LE[1], RE[1], LE[2], RE[2]
        res = self.select_grid(LE, RE, level)
        cdef int eterm[3] 
        eterm[0] = eterm[1] = eterm[2] = 0
        if res == 0:
            for i in range(8):
                mask[root.local_ind,i] = 0
            return
        # Now we visit all our children.  We subtract off sdds for the first
        # pass because we center it on the first cell.
        spos[0] = pos[0] - sdds[0]/2.0
        cdef int next_level, this_level
        # next_level: an int that says whether or not we can progress to children
        # this_level: an int that says whether or not we can select from this
        # level
        next_level = this_level = 1
        if level == self.max_level:
            next_level = 0
        if level < self.min_level or level > self.max_level:
            this_level = 0
        for i in range(2):
            spos[1] = pos[1] - sdds[1]/2.0
            for j in range(2):
                spos[2] = pos[2] - sdds[2]/2.0
                for k in range(2):
                    ii = ((k*2)+j)*2+i
                    ch = root.children[i][j][k]
                    if this_level == 1 and ch != NULL:
                        mask[root.local_ind, ii] = 0
                        self.recursively_select_octs(
                            ch, spos, sdds, mask, level + 1)
                    elif this_level == 1:
                        mask[root.local_ind, ii] = \
                            self.select_cell(spos, sdds, eterm)
                    spos[2] += sdds[2]
                spos[1] += sdds[1]
            spos[0] += sdds[0]

    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level) nogil:
        return 0
    
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3],
                         int eterm[3]) nogil:
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
        cdef int i, j, k, ind[3][2]
        child_mask = gobj.child_mask
        for i in range(3):
            ind[i][0] = 0
            ind[i][1] = gobj.ActiveDimensions[i]
            dds[i] = odds[i]
        cdef int count = 0
        cdef int check = 1
        # Check for the level bounds
        cdef np.int32_t level = gobj.Level
        if level < self.min_level or level > self.max_level:
            return count
        # We set this to 1 if we ignore child_mask
        cdef int this_level = 0
        if level == self.max_level:
            this_level = 1
        cdef int eterm[3]
        self.set_bounds(<np.float64_t *> left_edge.data,
                        <np.float64_t *> right_edge.data,
                        dds, ind, &check)
        with nogil:
            if check == 1:
                pos[0] = left_edge[0] + dds[0] * 0.5
                for i in range(ind[0][0], ind[0][1]):
                    eterm[0] = 0
                    pos[1] = left_edge[1] + dds[1] * 0.5
                    for j in range(ind[1][0], ind[1][1]):
                        eterm[1] = 0
                        pos[2] = left_edge[2] + dds[2] * 0.5
                        for k in range(ind[2][0], ind[2][1]):
                            eterm[2] = 0
                            if child_mask[i,j,k] == 1 or this_level == 1:
                                count += self.select_cell(pos, dds, eterm)
                            if eterm[2] == 1: break
                            pos[2] += dds[1]
                        if eterm[1] == 1: break
                        pos[1] += dds[1]
                    if eterm[0] == 1: break
                    pos[0] += dds[0]
            else:
                for i in range(ind[0][0], ind[0][1]):
                    for j in range(ind[1][0], ind[1][1]):
                        for k in range(ind[2][0], ind[2][1]):
                            if child_mask[i,j,k] == 1 or this_level == 1:
                                count += 1
        return count

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        child_mask = gobj.child_mask
        cdef np.ndarray[np.uint8_t, ndim=3] mask 
        cdef int ind[3][2]
        cdef np.ndarray[np.float64_t, ndim=1] odds = gobj.dds
        cdef np.ndarray[np.float64_t, ndim=1] left_edge = gobj.LeftEdge
        cdef np.ndarray[np.float64_t, ndim=1] right_edge = gobj.RightEdge
        cdef int i, j, k
        cdef np.float64_t dds[3], pos[3]
        for i in range(3):
            dds[i] = odds[i]
            ind[i][0] = 0
            ind[i][1] = gobj.ActiveDimensions[i]
        mask = np.zeros(gobj.ActiveDimensions, dtype='uint8')
        cdef int check = 1
        cdef int total = 0
        cdef int eterm[3]
        self.set_bounds(<np.float64_t *> left_edge.data,
                        <np.float64_t *> right_edge.data,
                        dds, ind, &check)
        cdef int temp
        # Check for the level bounds
        cdef np.int32_t level = gobj.Level
        if level < self.min_level or level > self.max_level:
            return mask
        # We set this to 1 if we ignore child_mask
        cdef int this_level = 0
        if level == self.max_level:
            this_level = 1
        with nogil:
            if check == 1:
                pos[0] = left_edge[0] + dds[0] * 0.5
                for i in range(ind[0][0], ind[0][1]):
                    eterm[0] = 0
                    pos[1] = left_edge[1] + dds[1] * 0.5
                    for j in range(ind[1][0], ind[1][1]):
                        eterm[1] = 0
                        pos[2] = left_edge[2] + dds[2] * 0.5
                        for k in range(ind[2][0], ind[2][1]):
                            eterm[2] = 0
                            if child_mask[i,j,k] == 1 or this_level == 1:
                                mask[i,j,k] = self.select_cell(pos, dds, eterm)
                                total += mask[i,j,k]
                            if eterm[2] == 1: break
                            pos[2] += dds[1]
                        if eterm[1] == 1: break
                        pos[1] += dds[1]
                    if eterm[0] == 1: break
                    pos[0] += dds[0]
            else:
                for i in range(ind[0][0], ind[0][1]):
                    for j in range(ind[1][0], ind[1][1]):
                        for k in range(ind[2][0], ind[2][1]):
                            if child_mask[i,j,k] == 1 or this_level == 1:
                                mask[i,j,k] = 1
                            total += mask[i,j,k]
        if total == 0: return None
        return mask.astype("bool")

    cdef void set_bounds(self,
                         np.float64_t left_edge[3], np.float64_t right_edge[3],
                         np.float64_t dds[3], int ind[3][2], int *check):
        return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_points(self, np.ndarray[np.float64_t, ndim=1] x,
                           np.ndarray[np.float64_t, ndim=1] y,
                           np.ndarray[np.float64_t, ndim=1] z,
                           np.float64_t radius = 0.0):
        cdef int count = 0
        cdef int i
        cdef np.float64_t dds[3], pos[3]
        dds[0] = dds[1] = dds[2] = radius
        cdef int eterm[3]
        with nogil:
            for i in range(x.shape[0]):
                pos[0] = x[i]
                pos[1] = y[i]
                pos[2] = z[i]
                count += self.select_cell(pos, dds, eterm)
        return count

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_points(self, np.ndarray[np.float64_t, ndim=1] x,
                            np.ndarray[np.float64_t, ndim=1] y,
                            np.ndarray[np.float64_t, ndim=1] z,
                            np.float64_t radius = 0.0):
        cdef int count = 0
        cdef int i
        cdef np.float64_t dds[3], pos[3]
        dds[0] = dds[1] = dds[2] = radius
        cdef int eterm[3]
        cdef np.ndarray[np.uint8_t, ndim=1] mask 
        mask = np.zeros(x.shape[0], dtype='uint8')
        with nogil:
            for i in range(x.shape[0]):
                pos[0] = x[i]
                pos[1] = y[i]
                pos[2] = z[i]
                mask[i] = self.select_cell(pos, dds, eterm)
                count += mask[i]
        if count == 0: return None
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
                               np.float64_t right_edge[3],
                               np.int32_t level) nogil:
        cdef np.float64_t box_center, relcenter, closest, dist, edge
        return 1
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
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3],
                         int eterm[3]) nogil:
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
                               np.float64_t right_edge[3],
                               np.int32_t level) nogil:
        if level < self.min_level or level > self.max_level: return 0
        for i in range(3):
            if left_edge[i] >= self.right_edge[i]: return 0
            if right_edge[i] <= self.left_edge[i]: return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3],
                         int eterm[3]) nogil:
        cdef np.float64_t dxp
        for i in range(3):
            dxp = self.dx_pad * dds[i]
            if pos[i] - dxp >= self.right_edge[i]:
                eterm[i] = 1
                return 0
            if pos[i] + dxp <= self.left_edge[i]: return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void set_bounds(self,
                         np.float64_t left_edge[3], np.float64_t right_edge[3],
                         np.float64_t dds[3], int ind[3][2], int *check):
        cdef int temp, i, all_inside
        # Left pos is left_edge + 0.5 * dds
        all_inside = 1
        for i in range(3):
            if self.left_edge[i] > left_edge[i]:
                temp = <int> ((self.left_edge[i] - left_edge[i])/dds[i]) - 1
                ind[i][0] = iclip(temp, ind[i][0], ind[i][1] - 1)
                all_inside = 0
            if self.right_edge[i] < right_edge[i]:
                temp = <int> ((self.right_edge[i] - right_edge[i])/dds[i]) + 1
                ind[i][1] = iclip(temp, ind[i][0] + 1, ind[i][1])
                all_inside = 0
        check[0] = all_inside


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
                               np.float64_t right_edge[3],
                               np.int32_t level) nogil:
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
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3],
                         int eterm[3]) nogil:
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

cdef class CuttingPlaneSelector(SelectorObject):
    cdef np.float64_t norm_vec[3]
    cdef np.float64_t d

    def __init__(self, dobj):
        cdef int i
        for i in range(3):
            self.norm_vec[i] = dobj._norm_vec[i]
        self.d = dobj._d

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level) nogil:
        cdef int i, j, k, n
        cdef np.float64_t *arr[2]
        cdef np.float64_t pos[3]
        cdef np.float64_t gd
        arr[0] = left_edge
        arr[1] = right_edge
        all_under = 1
        all_over = 1
        # Check each corner
        for i in range(2):
            pos[0] = arr[i][0]
            for j in range(2):
                pos[1] = arr[j][1]
                for k in range(2):
                    pos[2] = arr[k][2]
                    gd = self.d
                    for n in range(3):
                        gd += pos[n] * self.norm_vec[n]
                    if gd < 0: all_over = 0
                    if gd > 0: all_under = 0
        if all_over == 1 or all_under == 1:
            return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3],
                         int eterm[3]) nogil:
        cdef np.float64_t diag2, height
        cdef int i
        height = self.d
        diag2 = 0
        for i in range(3):
            height += pos[i] * self.norm_vec[i]
            diag2 += dds[i] * dds[i] * 0.25
        if height * height <= diag2: return 1
        return 0

cutting_selector = CuttingPlaneSelector

cdef class SliceSelector(SelectorObject):
    cdef int axis
    cdef np.float64_t coord

    def __init__(self, dobj):
        self.axis = dobj.axis
        self.coord = dobj.coord

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void set_bounds(self,
                         np.float64_t left_edge[3], np.float64_t right_edge[3],
                         np.float64_t dds[3], int ind[3][2], int *check):
        cdef int i
        for i in range(3):
            if self.axis == i:
                ind[i][0] = <int> ((self.coord - left_edge[i])/dds[i])
                ind[i][1] = ind[i][0] + 1
        check[0] = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level) nogil:
        if right_edge[self.axis] > self.coord \
           and left_edge[self.axis] <= self.coord:
            return 1
        return 0
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3],
                         int eterm[3]) nogil:
        if pos[self.axis] + 0.5*dds[self.axis] > self.coord \
           and pos[self.axis] - 0.5*dds[self.axis] <= self.coord:
            return 1
        return 0

slice_selector = SliceSelector

cdef class OrthoRaySelector(SelectorObject):

    cdef np.uint8_t px_ax
    cdef np.uint8_t py_ax
    cdef np.float64_t px
    cdef np.float64_t py
    cdef int axis

    def __init__(self, dobj):
        self.axis = dobj.axis
        self.px_ax = dobj.px_ax
        self.py_ax = dobj.py_ax
        self.px = dobj.px
        self.py = dobj.py

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level) nogil:
        if (    (self.px >= left_edge[self.px_ax])
            and (self.px < right_edge[self.px_ax])
            and (self.py >= left_edge[self.py_ax])
            and (self.py < right_edge[self.py_ax])):
            return 1
        return 0
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3],
                         int eterm[3]) nogil:
        if (    (self.px >= pos[self.px_ax] - 0.5*dds[self.px_ax])
            and (self.px <  pos[self.px_ax] + 0.5*dds[self.px_ax])
            and (self.py >= pos[self.py_ax] - 0.5*dds[self.py_ax])
            and (self.py <  pos[self.py_ax] + 0.5*dds[self.py_ax])):
            return 1
        return 0


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void set_bounds(self,
                         np.float64_t left_edge[3], np.float64_t right_edge[3],
                         np.float64_t dds[3], int ind[3][2], int *check):
        cdef int i
        for i in range(3):
            if self.px_ax == i:
                ind[i][0] = <int> ((self.px - left_edge[i])/dds[i])
                ind[i][1] = ind[i][0] + 1
            elif self.py_ax == i:
                ind[i][0] = <int> ((self.py - left_edge[i])/dds[i])
                ind[i][1] = ind[i][0] + 1
        check[0] = 0

ortho_ray_selector = OrthoRaySelector

cdef struct IntegrationAccumulator:
    np.float64_t *t
    np.float64_t *dt
    np.uint8_t *child_mask
    int hits

cdef void dt_sampler(
             VolumeContainer *vc,
             np.float64_t v_pos[3],
             np.float64_t v_dir[3],
             np.float64_t enter_t,
             np.float64_t exit_t,
             int index[3],
             void *data) nogil:
    cdef IntegrationAccumulator *am = <IntegrationAccumulator *> data
    cdef int di = (index[0]*vc.dims[1]+index[1])*vc.dims[2]+index[2] 
    if am.child_mask[di] == 0 or enter_t == exit_t:
        return
    am.hits += 1
    am.t[di] = enter_t
    am.dt[di] = (exit_t - enter_t)

cdef class RaySelector(SelectorObject):

    cdef np.float64_t p1[3]
    cdef np.float64_t p2[3]
    cdef np.float64_t vec[3]

    def __init__(self, dobj):
        cdef int i
        for i in range(3):
            self.vec[i] = dobj.vec[i]
            self.p1[i] = dobj.start_point[i]
            self.p2[i] = dobj.end_point[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level) nogil:
        cdef int i, ax
        cdef int i1, i2
        cdef np.float64_t vs[3], t, v[3]
        for ax in range(3):
            i1 = (ax+1) % 3
            i2 = (ax+2) % 3
            t = (left_edge[ax] - self.p1[ax])/self.vec[ax]
            for i in range(3):
                vs[i] = t * self.vec[i] + self.p1[i]
            if left_edge[i1] <= vs[i1] and \
               right_edge[i1] >= vs[i1] and \
               left_edge[i2] <= vs[i2] and \
               right_edge[i2] >= vs[i2] and \
               0.0 <= t <= 1.0:
                return 1
            t = (right_edge[ax] - self.p1[ax])/self.vec[ax]
            for i in range(3):
                vs[i] = t * self.vec[i] + self.p1[i]
            if left_edge[i1] <= vs[i1] and \
               right_edge[i1] >= vs[i1] and \
               left_edge[i2] <= vs[i2] and \
               right_edge[i2] >= vs[i2] and\
               0.0 <= t <= 1.0:
                return 1
        # if the point is fully enclosed, we count the grid
        if left_edge[0] <= self.p1[0] and \
           right_edge[0] >= self.p1[0] and \
           left_edge[1] <= self.p1[1] and \
           right_edge[1] >= self.p1[1] and \
           left_edge[2] <= self.p1[2] and \
           right_edge[2] >= self.p1[2]:
            return 1
        if left_edge[0] <= self.p2[0] and \
           right_edge[0] >= self.p2[0] and \
           left_edge[1] <= self.p2[1] and \
           right_edge[1] >= self.p2[1] and \
           left_edge[2] <= self.p2[2] and \
           right_edge[2] >= self.p2[2]:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.float64_t, ndim=3] t, dt
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef int i
        cdef IntegrationAccumulator ia
        cdef VolumeContainer vc
        mask = np.zeros(gobj.ActiveDimensions, dtype='uint8')
        t = np.zeros(gobj.ActiveDimensions, dtype="float64")
        dt = np.zeros(gobj.ActiveDimensions, dtype="float64") - 1
        child_mask = gobj.child_mask
        ia.t = <np.float64_t *> t.data
        ia.dt = <np.float64_t *> dt.data
        ia.child_mask = <np.uint8_t *> child_mask.data
        ia.hits = 0
        for i in range(3):
            vc.left_edge[i] = gobj.LeftEdge[i]
            vc.right_edge[i] = gobj.RightEdge[i]
            vc.dds[i] = gobj.dds[i]
            vc.idds[i] = 1.0/gobj.dds[i]
            vc.dims[i] = dt.shape[i]
        walk_volume(&vc, self.p1, self.vec, dt_sampler, <void*> &ia)
        for i in range(dt.shape[0]):
            for j in range(dt.shape[1]):
                for k in range(dt.shape[2]):
                    if dt[i,j,k] >= 0:
                        mask[i,j,k] = 1
        return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_cells(self, gobj):
        return self.fill_mask(gobj).sum()
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def get_dt(self, gobj):
        cdef np.ndarray[np.float64_t, ndim=3] t, dt
        cdef np.ndarray[np.float64_t, ndim=1] tr, dtr
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef int i, j, k, ni
        cdef IntegrationAccumulator ia
        cdef VolumeContainer vc
        t = np.zeros(gobj.ActiveDimensions, dtype="float64")
        dt = np.zeros(gobj.ActiveDimensions, dtype="float64") - 1
        child_mask = gobj.child_mask
        ia.t = <np.float64_t *> t.data
        ia.dt = <np.float64_t *> dt.data
        ia.child_mask = <np.uint8_t *> child_mask.data
        ia.hits = 0
        for i in range(3):
            vc.left_edge[i] = gobj.LeftEdge[i]
            vc.right_edge[i] = gobj.RightEdge[i]
            vc.dds[i] = gobj.dds[i]
            vc.idds[i] = 1.0/gobj.dds[i]
            vc.dims[i] = dt.shape[i]
        walk_volume(&vc, self.p1, self.vec, dt_sampler, <void*> &ia)
        tr = np.zeros(ia.hits, dtype="float64")
        dtr = np.zeros(ia.hits, dtype="float64")
        ni = 0
        for i in range(dt.shape[0]):
            for j in range(dt.shape[1]):
                for k in range(dt.shape[2]):
                    if dt[i,j,k] >= 0:
                        tr[ni] = t[i,j,k]
                        dtr[ni] = dt[i,j,k]
                        ni += 1
        if not (ni == ia.hits):
            print ni, ia.hits
        return dtr, tr
    
ray_selector = RaySelector

cdef class DataCollectionSelector(SelectorObject):
    cdef object obj_ids
    cdef np.int64_t nids

    def __init__(self, dobj):
        self.obj_ids = dobj._obj_ids
        self.nids = self.obj_ids.shape[0]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void set_bounds(self,
                         np.float64_t left_edge[3], np.float64_t right_edge[3],
                         np.float64_t dds[3], int ind[3][2], int *check):
        check[0] = 0
        return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        cdef int i, n
        cdef int ng = left_edges.shape[0]
        cdef np.ndarray[np.uint8_t, ndim=1] gridi = np.zeros(ng, dtype='uint8')
        cdef np.ndarray[np.int64_t, ndim=1] oids = self.obj_ids
        with nogil:
            for n in range(self.nids):
                # Call our selector function
                # Check if the sphere is inside the grid
                gridi[oids[n]] = 1
        return gridi.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_cells(self, gobj):
        return gobj.ActiveDimensions.prod()
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.uint8_t, ndim=3] mask 
        mask = np.ones(gobj.ActiveDimensions, dtype='uint8')
        return mask.astype("bool")

data_collection_selector = DataCollectionSelector

cdef class EllipsoidSelector(SelectorObject):
    cdef np.float64_t vec[3][3]
    cdef np.float64_t mag[3]
    cdef np.float64_t center[3]

    def __init__(self, dobj):
        cdef int i
        for i in range(3):
            self.center[i] = dobj.center[i]
            self.vec[0][i] = dobj._e0[i]
            self.vec[1][i] = dobj._e1[i]
            self.vec[2][i] = dobj._e2[i]
        self.mag[0] = dobj._A
        self.mag[1] = dobj._B
        self.mag[2] = dobj._C

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level) nogil:
        # This is the sphere selection
        cdef np.float64_t radius2, box_center, relcenter, closest, dist, edge
        return 1
        radius2 = self.mag[0] * self.mag[0]
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
        if dist < radius2: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3],
                         int eterm[3]) nogil:
        cdef np.float64_t dot_evec[3]
        cdef np.float64_t dist
        cdef int i, j
        dot_evec[0] = dot_evec[1] = dot_evec[2] = 0
        # Calculate the rotated dot product
        for i in range(3): # axis
            dist = pos[i] - self.center[i]
            for j in range(3):
                dot_evec[j] += dist * self.vec[j][i]
        dist = 0.0
        for i in range(3):
            dist += (dot_evec[i] * dot_evec[i])/(self.mag[i] * self.mag[i])
        if dist <= 1.0: return 1
        return 0

ellipsoid_selector = EllipsoidSelector

cdef class GridSelector(SelectorObject):
    cdef object ind

    def __init__(self, dobj):
        self.ind = dobj.id - dobj._id_offset

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void set_bounds(self,
                         np.float64_t left_edge[3], np.float64_t right_edge[3],
                         np.float64_t dds[3], int ind[3][2], int *check):
        check[0] = 0
        return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        cdef int ng = left_edges.shape[0]
        cdef np.ndarray[np.uint8_t, ndim=1] gridi = np.zeros(ng, dtype='uint8')
        gridi[self.ind] = 1
        return gridi.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_cells(self, gobj):
        return gobj.ActiveDimensions.prod()
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        return np.ones(gobj.ActiveDimensions, dtype='bool')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3],
                         int eterm[3]) nogil:
        return 1


grid_selector = GridSelector

