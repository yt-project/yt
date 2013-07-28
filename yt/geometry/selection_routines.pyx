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
from libc.stdlib cimport malloc, free
from fp_utils cimport fclip, iclip
from selection_routines cimport SelectorObject
from oct_container cimport OctreeContainer, OctAllocationContainer, Oct
cimport oct_visitors
from oct_visitors cimport cind
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

cdef class SelectorObject:

    def __cinit__(self, dobj):
        self.min_level = getattr(dobj, "min_level", 0)
        self.max_level = getattr(dobj, "max_level", 99)
        self.overlap_cells = 0

        for i in range(3) :
            if dobj.pf.periodicity[i] and dobj.pf.domain_left_edge[i] != 0.0 :
                print "SelectorObject periodicity assumes left_edge == 0"
                raise RuntimeError

            self.domain_width[i] = dobj.pf.domain_right_edge[i] - \
                                   dobj.pf.domain_left_edge[i]
            self.periodicity[i] = dobj.pf.periodicity[i]

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

    def count_octs(self, OctreeContainer octree, int domain_id = -1):
        cdef OctVisitorData data
        data.index = 0
        data.last = -1
        data.domain = domain_id
        octree.visit_all_octs(self, oct_visitors.count_total_octs, &data)
        return data.index

    def count_oct_cells(self, OctreeContainer octree, int domain_id = -1):
        cdef OctVisitorData data
        data.index = 0
        data.domain = domain_id
        octree.visit_all_octs(self, oct_visitors.count_total_cells, &data)
        return data.index

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void recursively_visit_octs(self, Oct *root,
                        np.float64_t pos[3], np.float64_t dds[3],
                        int level, 
                        oct_visitor_function *func,
                        OctVisitorData *data,
                        int visit_covered = 0):
        # visit_covered tells us whether this octree supports partial
        # refinement.  If it does, we need to handle this specially -- first
        # we visit *this* oct, then we make a second pass to check any child
        # octs.
        cdef np.float64_t LE[3], RE[3], sdds[3], spos[3]
        cdef int i, j, k, res, mi
        cdef Oct *ch
        cdef np.uint8_t selected
        # Remember that pos is the *center* of the oct, and dds is the oct
        # width.  So to get to the edges, we add/subtract half of dds.
        for i in range(3):
            # sdds is the cell width
            sdds[i] = dds[i]/2.0
            LE[i] = pos[i] - dds[i]/2.0
            RE[i] = pos[i] + dds[i]/2.0
        #print LE[0], RE[0], LE[1], RE[1], LE[2], RE[2]
        res = self.select_grid(LE, RE, level, root)
        if res == 1 and data.domain > 0 and root.domain != data.domain:
            res = -1
        cdef int increment = 1
        cdef int next_level, this_level
        # next_level: an int that says whether or not we can progress to children
        # this_level: an int that says whether or not we can select from this
        # level
        next_level = this_level = 1
        if res == -1:
            # This happens when we do domain selection but the oct has
            # children.  This would allow an oct to pass to its children but
            # not get accessed itself.
            next_level = 1
            this_level = 0
        elif level == self.max_level:
            next_level = 0
        elif level < self.min_level or level > self.max_level:
            this_level = 0
        if res == 0 and this_level == 1:
            return
        # Now we visit all our children.  We subtract off sdds for the first
        # pass because we center it on the first cell.
        cdef int iter = 1 - visit_covered # 2 if 1, 1 if 0.
        while iter < 2:
            spos[0] = pos[0] - sdds[0]/2.0
            for i in range(2):
                spos[1] = pos[1] - sdds[1]/2.0
                for j in range(2):
                    spos[2] = pos[2] - sdds[2]/2.0
                    for k in range(2):
                        ch = NULL
                        if root.children != NULL:
                            ch = root.children[cind(i,j,k)]
                        if iter == 1 and next_level == 1 and ch != NULL:
                            data.pos[0] = (data.pos[0] << 1) + i
                            data.pos[1] = (data.pos[1] << 1) + j
                            data.pos[2] = (data.pos[2] << 1) + k
                            data.level += 1
                            self.recursively_visit_octs(
                                ch, spos, sdds, level + 1, func, data,
                                visit_covered)
                            data.pos[0] = (data.pos[0] >> 1)
                            data.pos[1] = (data.pos[1] >> 1)
                            data.pos[2] = (data.pos[2] >> 1)
                            data.level -= 1
                        elif this_level == 1:
                            selected = self.select_cell(spos, sdds)
                            if ch != NULL:
                                selected *= self.overlap_cells
                            data.global_index += increment
                            increment = 0
                            data.ind[0] = i
                            data.ind[1] = j
                            data.ind[2] = k
                            func(root, data, selected)
                        spos[2] += sdds[2]
                    spos[1] += sdds[1]
                spos[0] += sdds[0]
            this_level = 0 # We turn this off for the second pass.
            iter += 1

    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level, Oct *o = NULL) nogil:
        return 0
    
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 0

    cdef int select_point(self, np.float64_t pos[3] ) nogil:
        return 0

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius ) nogil:
        return 0

    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.float64_t difference(self, np.float64_t x1, np.float64_t x2, int d) nogil:
        cdef np.float64_t rel = x1 - x2
        if self.periodicity[d] :
            if rel > self.domain_width[d]/2.0 :
                rel -= self.domain_width[d]
            elif rel < -self.domain_width[d]/2.0 :
                rel += self.domain_width[d]
        return rel

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_cells(self, gobj):
        cdef np.ndarray[np.float64_t, ndim=1] odds = gobj.dds
        cdef np.ndarray[np.float64_t, ndim=1] left_edge = gobj.LeftEdge
        cdef np.ndarray[np.float64_t, ndim=1] right_edge = gobj.RightEdge
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef np.float64_t dds[3], pos[3]
        cdef int i, j, k, dim[3]
        child_mask = gobj.child_mask
        for i in range(3):
            dim[i] = gobj.ActiveDimensions[i]
            dds[i] = odds[i]
        cdef int count = 0
        # Check for the level bounds
        cdef np.int32_t level = gobj.Level
        if level < self.min_level or level > self.max_level:
            return count
        # We set this to 1 if we ignore child_mask
        cdef int this_level = 0
        if level == self.max_level:
            this_level = 1
        with nogil:
            pos[0] = left_edge[0] + dds[0] * 0.5
            for i in range(dim[0]):
                pos[1] = left_edge[1] + dds[1] * 0.5
                for j in range(dim[1]):
                    pos[2] = left_edge[2] + dds[2] * 0.5
                    for k in range(dim[2]):
                        if child_mask[i,j,k] == 1 or this_level == 1:
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
        cdef int ind[3][2]
        cdef int dim[3]
        cdef np.ndarray[np.float64_t, ndim=1] odds = gobj.dds
        cdef np.ndarray[np.float64_t, ndim=1] left_edge = gobj.LeftEdge
        cdef np.ndarray[np.float64_t, ndim=1] right_edge = gobj.RightEdge
        cdef int i, j, k
        cdef np.float64_t dds[3], pos[3]
        for i in range(3):
            dds[i] = odds[i]
            dim[i] = gobj.ActiveDimensions[i]
        mask = np.zeros(gobj.ActiveDimensions, dtype='uint8')
        cdef int total = 0
        cdef int temp
        # Check for the level bounds
        cdef np.int32_t level = gobj.Level
        if level < self.min_level or level > self.max_level:
            return mask.astype("bool")
        # We set this to 1 if we ignore child_mask
        cdef int this_level = 0
        if level == self.max_level:
            this_level = 1
        with nogil:
            pos[0] = left_edge[0] + dds[0] * 0.5
            for i in range(dim[0]):
                pos[1] = left_edge[1] + dds[1] * 0.5
                for j in range(dim[1]):
                    pos[2] = left_edge[2] + dds[2] * 0.5
                    for k in range(dim[2]):
                        if child_mask[i,j,k] == 1 or this_level == 1:
                            mask[i,j,k] = self.select_cell(pos, dds)
                            total += mask[i,j,k]
                        pos[2] += dds[1]
                    pos[1] += dds[1]
                pos[0] += dds[0]
        if total == 0: return None
        return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_points(self, np.ndarray[np.float64_t, ndim=1] x,
                           np.ndarray[np.float64_t, ndim=1] y,
                           np.ndarray[np.float64_t, ndim=1] z,
                           np.float64_t radius = 0.0):
        cdef int count = 0
        cdef int i
        cdef np.float64_t pos[3]
        with nogil:
            if radius == 0.0 :
                for i in range(x.shape[0]):
                    pos[0] = x[i]
                    pos[1] = y[i]
                    pos[2] = z[i]
                    count += self.select_point(pos)
            else :
                for i in range(x.shape[0]):
                    pos[0] = x[i]
                    pos[1] = y[i]
                    pos[2] = z[i]
                    count += self.select_sphere(pos, radius)
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
        cdef np.float64_t pos[3]
        cdef np.ndarray[np.uint8_t, ndim=1] mask 
        mask = np.zeros(x.shape[0], dtype='uint8')

        # this is to allow selectors to optimize the point vs
        # 0-radius sphere case.  These two may have different 
        # effects for 0-volume selectors, however (collision 
        # between a ray and a point is null, while ray and a
        # sphere is allowed)
        with nogil:
            if radius == 0.0 :
                for i in range(x.shape[0]) :
                    pos[0] = x[i]
                    pos[1] = y[i]
                    pos[2] = z[i]
                    mask[i] = self.select_point(pos)
                    count += mask[i]
            else :
                for i in range(x.shape[0]):
                    pos[0] = x[i]
                    pos[1] = y[i]
                    pos[2] = z[i]
                    mask[i] = self.select_sphere(pos, radius)
                    count += mask[i]
        if count == 0: return None
        return mask.astype("bool")

    def __hash__(self):
        return hash(self._hash_vals() + self._base_hash())

    def _hash_vals(self):
        raise NotImplementedError

    def _base_hash(self):
        return (self.min_level, self.max_level, self.overlap_cells,
                self.periodicity[0],
                self.periodicity[1],
                self.periodicity[2],
                self.domain_width[0],
                self.domain_width[1],
                self.domain_width[2])

cdef class SphereSelector(SelectorObject):
    cdef np.float64_t radius
    cdef np.float64_t radius2
    cdef np.float64_t center[3]

    def __init__(self, dobj):
        for i in range(3):
            self.center[i] = dobj.center[i]
        self.radius = dobj.radius
        self.radius2 = dobj.radius * dobj.radius

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level, Oct *o = NULL) nogil:
        return self.select_bbox(left_edge,right_edge)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        # sphere center either inside cell or center of cell lies inside sphere
        if (pos[0] - 0.5*dds[0] <= self.center[0] <= pos[0]+0.5*dds[0] and
            pos[1] - 0.5*dds[1] <= self.center[1] <= pos[1]+0.5*dds[1] and
            pos[2] - 0.5*dds[2] <= self.center[2] <= pos[2]+0.5*dds[2]):
            return 1
        return self.select_point(pos)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int i
        cdef np.float64_t dist, dist2 = 0
        for i in range(3):
            dist = self.difference( pos[i], self.center[i], i )
            dist2 += dist*dist
        if dist2 <= self.radius2: return 1
        return 0
   
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int i
        cdef np.float64_t dist, dist2 = 0
        for i in range(3):
            dist = self.difference( pos[i], self.center[i], i ) 
            dist2 += dist*dist
        dist = self.radius+radius
        if dist2 <= dist*dist: return 1
        return 0
 
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef np.float64_t box_center, relcenter, closest, dist, edge
        cdef int i
        if (left_edge[0] <= self.center[0] <= right_edge[0] and
            left_edge[1] <= self.center[1] <= right_edge[1] and
            left_edge[2] <= self.center[2] <= right_edge[2]):
            return 1
        # http://www.gamedev.net/topic/335465-is-this-the-simplest-sphere-aabb-collision-test/
        dist = 0
        for i in range(3):
            box_center = (right_edge[i] + left_edge[i])/2.0
            relcenter = self.difference( box_center, self.center[i], i )
            edge = right_edge[i] - left_edge[i]
            closest = relcenter - fclip(relcenter, -edge/2.0, edge/2.0)
            dist += closest*closest
        if dist <= self.radius2: return 1
        return 0

    def _hash_vals(self):
        return (self.radius, self.radius2,
                self.center[0], self.center[1], self.center[2])

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
                               np.int32_t level, Oct *o = NULL) nogil:
        if level < self.min_level or level > self.max_level: return 0
        return self.select_bbox( left_edge, right_edge )

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int i, shift, included
        cdef np.float64_t LE, RE
        for i in range(3):
            if not self.periodicity[i]:
                if left_edge[i] > self.right_edge[i]: return 0
                if right_edge[i] < self.left_edge[i]: return 0
            else:
                included = 1
                for shift in range(3):
                    LE = left_edge[i] + self.domain_width[i] * (shift - 1)
                    RE = right_edge[i] + self.domain_width[i] * (shift - 1)
                    if LE > self.right_edge[i]: continue
                    if RE < self.left_edge[i]: continue
                    included = 1
                if included == 0: return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef int i
        cdef np.float64_t LE[3], RE[3]
        if self.dx_pad == 0.0:
            return self.select_point(pos)
        for i in range(3):
            LE[i] = pos[i] - self.dx_pad * dds[i]
            RE[i] = pos[i] + self.dx_pad * dds[i]
        return self.select_bbox(LE, RE)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        # assume pos[3] is inside domain
        cdef int i, shift, included
        cdef np.float64_t dxp, ppos
        for i in range(3):
            if not self.periodicity[i]:
                if pos[i] < self.left_edge[i]: return 0
                if pos[i] > self.right_edge[i]: return 0
            else:
                included = 0
                for shift in range(3):
                    ppos = pos[i] + self.domain_width[i] * (shift - 1)
                    if ppos < self.left_edge[i]: continue
                    if ppos > self.right_edge[i]: continue
                    included = 1
                if included == 0: return 0
        return 1

    def _hash_vals(self):
        return (self.left_edge[0], self.left_edge[1], self.left_edge[2],
                self.right_edge[0], self.right_edge[1], self.right_edge[2],
                self.dx_pad)

region_selector = RegionSelector

cdef class DiskSelector(SelectorObject):
    cdef np.float64_t norm_vec[3]
    cdef np.float64_t center[3]
    cdef np.float64_t radius, radius2
    cdef np.float64_t height

    def __init__(self, dobj):
        cdef int i
        for i in range(3):
            self.norm_vec[i] = dobj._norm_vec[i]
            self.center[i] = dobj.center[i]
        self.radius = dobj._radius
        self.radius2 = dobj._radius * dobj._radius
        self.height = dobj._height

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level, Oct *o = NULL) nogil:
        return self.select_bbox( left_edge, right_edge )

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return self.select_point( pos ) 

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef np.float64_t h, d, r2, temp
        cdef int i
        h = d = 0
        for i in range(3):
            temp = self.difference( pos[i], self.center[i], i )
            h += temp * self.norm_vec[i]
            d += temp*temp
        r2 = (d - h*h)
        if fabs(h) <= self.height and r2 <= self.radius2: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius ) nogil:
        cdef np.float64_t h, d, r2, temp
        cdef int i
        h = d = 0
        for i in range(3):
            temp = self.difference( pos[i], self.center[i], i )
            h += pos[i] * self.norm_vec[i]
            d += temp*temp
        r2 = (d - h*h)
        d = self.radius+radius
        if fabs(h) <= self.height+radius and r2 <= d*d: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3] ) nogil:
        cdef np.float64_t *arr[2]
        cdef np.float64_t pos[3], H, D, R2, temp
        cdef int i, j, k, n
        cdef int all_under = 1
        cdef int all_over = 1
        cdef int any_radius = 0
        # A moment of explanation (revised):
        #    The disk and bounding box collide if any of the following are true:
        #    1) the center of the disk is inside the bounding box
        #    2) any corner of the box lies inside the disk
        #    3) the box spans the plane (!all_under and !all_over) and at least
        #       one corner is within the cylindrical radius

        # check if disk center lies inside bbox
        if left_edge[0] <= self.center[0] <= right_edge[0] and \
           left_edge[1] <= self.center[1] <= right_edge[1] and \
           left_edge[2] <= self.center[2] <= right_edge[2] :
            return 1
        
        # check all corners
        arr[0] = left_edge
        arr[1] = right_edge
        for i in range(2):
            pos[0] = arr[i][0]
            for j in range(2):
                pos[1] = arr[j][1]
                for k in range(2):
                    pos[2] = arr[k][2]
                    H = D = 0
                    for n in range(3):
                        temp = self.difference( pos[n], self.center[n], n )
                        H += (temp * self.norm_vec[n])
                        D += temp*temp
                    R2 = (D - H*H)
                    if R2 < self.radius2 :
                        any_radius = 1
                        if fabs(H) < self.height: return 1
                    if H < 0: all_over = 0
                    if H > 0: all_under = 0
        if all_over == 0 and all_under == 0 and any_radius == 1: return 1
        return 0

    def _hash_vals(self):
        return (self.norm_vec[0], self.norm_vec[1], self.norm_vec[2],
                self.center[0], self.center[1], self.center[2],
                self.radius, self.radius2, self.height)

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
                               np.int32_t level, Oct *o = NULL) nogil:
        return self.select_bbox(left_edge,right_edge)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef np.float64_t left_edge[3]
        cdef np.float64_t right_edge[3]
        cdef int i
        for i in range(3):
            left_edge[i] = pos[i] - 0.5*dds[i]
            right_edge[i] = pos[i] + 0.5*dds[i]
        return self.select_bbox(left_edge, right_edge)

    cdef int select_point(self, np.float64_t pos[3] ) nogil:
        # two 0-volume constructs don't intersect
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius ) nogil:
        cdef int i
        cdef np.float64_t height = self.d
        for i in range(3) :
            height += pos[i] * self.norm_vec[i]
        if height*height <= radius*radius : return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3] ) nogil:
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
                    # this allows corners and faces on the low-end to
                    # collide, while not selecting cells on the high-side 
                    if i == 0 and j == 0 and k == 0 :
                        if gd <= 0: all_over = 0
                        if gd >= 0: all_under = 0
                    else :
                        if gd < 0: all_over = 0
                        if gd > 0: all_under = 0
        if all_over == 1 or all_under == 1:
            return 0
        return 1

    def _hash_vals(self):
        return (self.norm_vec[0], self.norm_vec[1], self.norm_vec[2],
                self.d)

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
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level, Oct *o = NULL) nogil:
        return self.select_bbox( left_edge, right_edge )

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        if pos[self.axis] + 0.5*dds[self.axis] > self.coord \
           and pos[self.axis] - 0.5*dds[self.axis] <= self.coord:
            return 1
        return 0

    cdef int select_point(self, np.float64_t pos[3] ) nogil:
        # two 0-volume constructs don't intersect
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius ) nogil:
        cdef np.float64_t dist = self.difference( pos[self.axis], self.coord, self.axis )
        if dist*dist < radius*radius:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3] ) nogil:
        if left_edge[self.axis] <= self.coord < right_edge[self.axis]:
            return 1
        return 0

    def _hash_vals(self):
        return (self.axis, self.coord)

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
                               np.int32_t level, Oct *o = NULL) nogil:
        return self.select_bbox(left_edge,right_edge)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        if self.px >= pos[self.px_ax] - 0.5*dds[self.px_ax] and \
           self.px <  pos[self.px_ax] + 0.5*dds[self.px_ax] and \
           self.py >= pos[self.py_ax] - 0.5*dds[self.py_ax] and \
           self.py <  pos[self.py_ax] + 0.5*dds[self.py_ax]:
            return 1
        return 0

    cdef int select_point(self, np.float64_t pos[3] ) nogil:
        # two 0-volume constructs don't intersect
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius ) nogil:
        cdef np.float64_t dx = self.difference( pos[self.px_ax], self.px, self.px_ax )
        cdef np.float64_t dy = self.difference( pos[self.py_ax], self.py, self.py_ax )
        if dx*dx + dy*dy < radius*radius:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3] ) nogil:
        if left_edge[self.px_ax] <= self.px < right_edge[self.px_ax] and \
           left_edge[self.py_ax] <= self.py < right_edge[self.py_ax] :
            return 1
        return 0

    def _hash_vals(self):
        return (self.px_ax, self.py_ax, self.px, self.py, self.axis)

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
                               np.int32_t level, Oct *o = NULL) nogil:
        return self.select_bbox(left_edge,right_edge)

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

    cdef int select_point(self, np.float64_t pos[3] ) nogil:
        # two 0-volume constructs don't intersect
        return 0

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius ) nogil:
        # not implemented
        return 0        

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int i, ax
        cdef int i1, i2
        cdef np.float64_t vs[3], t, v[3]

        # if either point is fully enclosed, we select the bounding box
        if left_edge[0] <= self.p1[0] <= right_edge[0] and \
           left_edge[1] <= self.p1[1] <= right_edge[1] and \
           left_edge[2] <= self.p1[2] <= right_edge[2]:
            return 1
        if left_edge[0] <= self.p2[0] <= right_edge[0] and \
           left_edge[1] <= self.p2[1] <= right_edge[1] and \
           left_edge[2] <= self.p2[2] <= right_edge[2]:
            return 1

        for ax in range(3):
            i1 = (ax+1) % 3
            i2 = (ax+2) % 3
            t = (left_edge[ax] - self.p1[ax])/self.vec[ax]
            if 0.0 <= t <= 1.0 :
                for i in range(3):
                    vs[i] = t * self.vec[i] + self.p1[i]
                if left_edge[i1] <= vs[i1] <= right_edge[i1] and \
                   left_edge[i2] <= vs[i2] <= right_edge[i2] :
                    return 1
            t = (right_edge[ax] - self.p1[ax])/self.vec[ax]
            if 0.0 <= t <= 1.0 :
                for i in range(3):
                    vs[i] = t * self.vec[i] + self.p1[i]
                if left_edge[i1] <= vs[i1] <= right_edge[i1] and \
                   left_edge[i2] <= vs[i2] <= right_edge[i2] :
                    return 1
        return 0

    def _hash_vals(self):
        return (self.p1[0], self.p1[1], self.p1[2],
                self.p2[0], self.p2[1], self.p2[2],
                self.vec[0], self.vec[1], self.vec[2])

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

    def _hash_vals(self):
        return (hash(self.obj_ids.tostring()), self.nids)

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
                               np.int32_t level, Oct *o = NULL) nogil:
        return self.select_bbox(left_edge, right_edge)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return self.select_point(pos)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef np.float64_t dot_evec[3]
        cdef np.float64_t dist
        cdef int i, j
        dot_evec[0] = dot_evec[1] = dot_evec[2] = 0
        # Calculate the rotated dot product
        for i in range(3): # axis
            dist = self.difference(pos[i], self.center[i], i)
            for j in range(3):
                dot_evec[j] += dist * self.vec[j][i]
        dist = 0.0
        for i in range(3):
            dist += (dot_evec[i] * dot_evec[i])/(self.mag[i] * self.mag[i])
        if dist <= 1.0: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        # this is the sphere selection
        cdef int i
        cdef np.float64_t dist2 = 0
        for i in range(3):
            dist2 += self.difference( pos[i], self.center[i], i )**2
        if dist2 <= (self.mag[0]+radius)**2: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # This is the sphere selection
        cdef int i
        cdef np.float64_t box_center, relcenter, closest, dist, edge
        if left_edge[0] <= self.center[0] <= right_edge[0] and \
           left_edge[1] <= self.center[1] <= right_edge[1] and \
           left_edge[2] <= self.center[2] <= right_edge[2]:
            return 1
        # http://www.gamedev.net/topic/335465-is-this-the-simplest-sphere-aabb-collision-test/
        dist = 0
        for i in range(3):
            box_center = (right_edge[i] + left_edge[i])/2.0
            relcenter = self.difference(box_center, self.center[i], i)
            edge = right_edge[i] - left_edge[i]
            closest = relcenter - fclip(relcenter, -edge/2.0, edge/2.0)
            dist += closest * closest
        if dist <= self.mag[0]**2: return 1
        return 0

    def _hash_vals(self):
        return (self.vec[0][0], self.vec[0][1], self.vec[0][2],
                self.vec[1][0], self.vec[1][1], self.vec[1][2],
                self.vec[2][0], self.vec[2][1], self.vec[2][2],
                self.mag[0], self.mag[1], self.mag[2],
                self.center[0], self.center[1], self.center[2])

ellipsoid_selector = EllipsoidSelector

cdef class GridSelector(SelectorObject):
    cdef object ind

    def __init__(self, dobj):
        self.ind = dobj.id - dobj._id_offset

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
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 1

    cdef int select_point(self, np.float64_t pos[3] ) nogil:
        # we apparently don't check if the point actually lies in the grid..
        return 1

    def _hash_vals(self):
        return (self.ind,)

grid_selector = GridSelector

cdef class OctreeSubsetSelector(SelectorObject):
    cdef SelectorObject base_selector
    cdef public np.int64_t domain_id

    def __init__(self, dobj):
        self.base_selector = dobj.base_selector
        self.min_level = self.base_selector.min_level
        self.max_level = self.base_selector.max_level
        self.domain_id = dobj.domain_id
        self.overlap_cells = 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        raise RuntimeError

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        return self.base_selector.select_bbox(left_edge, right_edge)
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        # Because visitors now use select_grid, we should be explicitly
        # checking this.
        cdef int res
        res = self.base_selector.select_grid(left_edge, right_edge, level, o)
        if res == 1 and o != NULL and o.domain != self.domain_id:
            return -1
        return res
    
    def _hash_vals(self):
        return (hash(self.base_selector), self.domain_id)

octree_subset_selector = OctreeSubsetSelector

cdef class ParticleOctreeSubsetSelector(SelectorObject):
    # This is a numpy array, which will be a bool of ndim 1
    cdef np.uint64_t min_ind
    cdef np.uint64_t max_ind
    cdef SelectorObject base_selector

    def __init__(self, dobj):
        self.min_ind = dobj.min_ind
        self.max_ind = dobj.max_ind
        self.base_selector = dobj.base_selector
        self.min_level = self.base_selector.min_level
        self.max_level = self.base_selector.max_level

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        raise RuntimeError

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        return self.base_selector.select_bbox(left_edge, right_edge)
    
    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        # Because visitors now use select_grid, we should be explicitly
        # checking this.
        return self.base_selector.select_grid(left_edge, right_edge, level, o)
    
    def _hash_vals(self):
        return (hash(self.base_selector), self.min_ind, self.max_ind)

particle_octree_subset_selector = ParticleOctreeSubsetSelector

cdef class AlwaysSelector(SelectorObject):

    def __init__(self, dobj):
        self.overlap_cells = 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        cdef int ng = left_edges.shape[0]
        cdef np.ndarray[np.uint8_t, ndim=1] gridi = np.ones(ng, dtype='uint8')
        return gridi.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 1

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        return 1

    cdef int select_point(self, np.float64_t pos[3] ) nogil:
        return 1

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius ) nogil:
        return 1

    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        return 1

    def _hash_vals(self):
        return ("always", 1,)

always_selector = AlwaysSelector
