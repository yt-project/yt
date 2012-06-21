"""
AMR hierarchy container class

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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

import numpy as na

from yt.funcs import *
from yt.utilities.lib import \
    get_box_grids_level, \
    get_box_grids_below_level

class ObjectFindingMixin(object):

    def find_ray_grids(self, coord, axis):
        """
        Returns the (objects, indices) of grids that an (x,y) ray intersects
        along *axis*
        """
        # Let's figure out which grids are on the slice
        mask=na.ones(self.num_grids)
        # So if gRE > coord, we get a mask, if not, we get a zero
        #    if gLE > coord, we get a zero, if not, mask
        # Thus, if the coordinate is between the two edges, we win!
        na.choose(na.greater(self.grid_right_edge[:,x_dict[axis]],coord[0]),(0,mask),mask)
        na.choose(na.greater(self.grid_left_edge[:,x_dict[axis]],coord[0]),(mask,0),mask)
        na.choose(na.greater(self.grid_right_edge[:,y_dict[axis]],coord[1]),(0,mask),mask)
        na.choose(na.greater(self.grid_left_edge[:,y_dict[axis]],coord[1]),(mask,0),mask)
        ind = na.where(mask == 1)
        return self.grids[ind], ind

    def find_max(self, field, finest_levels = 3):
        """
        Returns (value, center) of location of maximum for a given field.
        """
        if (field, finest_levels) in self._max_locations:
            return self._max_locations[(field, finest_levels)]
        mg, mc, mv, pos = self.find_max_cell_location(field, finest_levels)
        self._max_locations[(field, finest_levels)] = (mv, pos)
        return mv, pos

    def find_max_cell_location(self, field, finest_levels = 3):
        if finest_levels is not False:
            gi = (self.grid_levels >= self.max_level - finest_levels).ravel()
            source = self.grid_collection([0.0]*3, self.grids[gi])
        else:
            source = self.all_data()
        mylog.debug("Searching %s grids for maximum value of %s",
                    len(source._grids), field)
        max_val, maxi, mx, my, mz, mg = \
            source.quantities["MaxLocation"]( field, lazy_reader=True)
        max_grid = self.grids[mg]
        mc = na.unravel_index(maxi, max_grid.ActiveDimensions)
        mylog.info("Max Value is %0.5e at %0.16f %0.16f %0.16f in grid %s at level %s %s", \
              max_val, mx, my, mz, max_grid, max_grid.Level, mc)
        self.parameters["Max%sValue" % (field)] = max_val
        self.parameters["Max%sPos" % (field)] = "%s" % ((mx,my,mz),)
        return max_grid, mc, max_val, na.array((mx,my,mz), dtype='float64')

    def find_min(self, field):
        """
        Returns (value, center) of location of minimum for a given field
        """
        gI = na.where(self.grid_levels >= 0) # Slow but pedantic
        minVal = 1e100
        for grid in self.grids[gI[0]]:
            mylog.debug("Checking %s (level %s)", grid.id, grid.Level)
            val, coord = grid.find_min(field)
            if val < minVal:
                minCoord = coord
                minVal = val
                minGrid = grid
        mc = na.array(minCoord)
        pos=minGrid.get_position(mc)
        mylog.info("Min Value is %0.5e at %0.16f %0.16f %0.16f in grid %s at level %s", \
              minVal, pos[0], pos[1], pos[2], minGrid, minGrid.Level)
        self.center = pos
        self.parameters["Min%sValue" % (field)] = minVal
        self.parameters["Min%sPos" % (field)] = "%s" % (pos)
        return minVal, pos

    def find_point(self, coord):
        """
        Returns the (objects, indices) of grids containing an (x,y,z) point
        """
        mask=na.ones(self.num_grids)
        for i in xrange(len(coord)):
            na.choose(na.greater(self.grid_left_edge[:,i],coord[i]), (mask,0), mask)
            na.choose(na.greater(self.grid_right_edge[:,i],coord[i]), (0,mask), mask)
        ind = na.where(mask == 1)
        return self.grids[ind], ind

    def find_field_value_at_point(self, fields, coord):
        r"""Find the value of fields at a point.
        
        Returns the values [field1, field2,...] of the fields at the given
        (x,y,z) point. Returns a list of field values in the same order
        as the input *fields*.
        
        Parameters
        ----------
        fields : string or list of strings
            The field(s) that will be returned.
        
        coord : list or array of floats
            The location for which field values will be returned.
        
        Examples
        --------
        >>> pf.h.find_field_value_at_point(['Density', 'Temperature'],
            [0.4, 0.3, 0.8])
        [2.1489e-24, 1.23843e4]
        """
        # Get the most-refined grid at this coordinate.
        this = self.find_point(coord)[0][-1]
        cellwidth = (this.RightEdge - this.LeftEdge) / this.ActiveDimensions
        mark = na.zeros(3).astype('int')
        # Find the index for the cell containing this point.
        for dim in xrange(len(coord)):
            mark[dim] = int((coord[dim] - this.LeftEdge[dim]) / cellwidth[dim])
        out = []
        fields = ensure_list(fields)
        # Pull out the values and add it to the out list.
        for field in fields:
            out.append(this[field][mark[0], mark[1], mark[2]])
        return out

    def find_slice_grids(self, coord, axis):
        """
        Returns the (objects, indices) of grids that a slice intersects along
        *axis*
        """
        # Let's figure out which grids are on the slice
        mask=na.ones(self.num_grids)
        # So if gRE > coord, we get a mask, if not, we get a zero
        #    if gLE > coord, we get a zero, if not, mask
        # Thus, if the coordinate is between the edges, we win!
        #ind = na.where( na.logical_and(self.grid_right_edge[:,axis] > coord, \
                                       #self.grid_left_edge[:,axis] < coord))
        na.choose(na.greater(self.grid_right_edge[:,axis],coord),(0,mask),mask)
        na.choose(na.greater(self.grid_left_edge[:,axis],coord),(mask,0),mask)
        ind = na.where(mask == 1)
        return self.grids[ind], ind

    def find_sphere_grids(self, center, radius):
        """
        Returns objects, indices of grids within a sphere
        """
        centers = (self.grid_right_edge + self.grid_left_edge)/2.0
        long_axis = na.maximum.reduce(self.grid_right_edge - self.grid_left_edge, 1)
        t = na.abs(centers - center)
        DW = self.parameter_file.domain_right_edge \
           - self.parameter_file.domain_left_edge
        na.minimum(t, na.abs(DW-t), t)
        dist = na.sqrt(na.sum((t**2.0), axis=1))
        gridI = na.where(dist < (radius + long_axis))
        return self.grids[gridI], gridI

    def get_box_grids(self, left_edge, right_edge):
        """
        Gets back all the grids between a left edge and right edge
        """
        grid_i = na.where((na.all(self.grid_right_edge > left_edge, axis=1)
                         & na.all(self.grid_left_edge < right_edge, axis=1)) == True)
        return self.grids[grid_i], grid_i

    def get_periodic_box_grids(self, left_edge, right_edge):
        mask = na.zeros(self.grids.shape, dtype='bool')
        dl = self.parameter_file.domain_left_edge
        dr = self.parameter_file.domain_right_edge
        left_edge = na.array(left_edge)
        right_edge = na.array(right_edge)
        dw = dr - dl
        left_dist = left_edge - dl
        db = right_edge - left_edge
        for off_x in [-1, 0, 1]:
            nle = left_edge.copy()
            nle[0] = (dw[0]*off_x + dl[0]) + left_dist[0]
            for off_y in [-1, 0, 1]:
                nle[1] = (dw[1]*off_y + dl[1]) + left_dist[1]
                for off_z in [-1, 0, 1]:
                    nle[2] = (dw[2]*off_z + dl[2]) + left_dist[2]
                    nre = nle + db
                    g, gi = self.get_box_grids(nle, nre)
                    mask[gi] = True
        return self.grids[mask], na.where(mask)

    def get_box_grids_below_level(self, left_edge, right_edge, level,
                                  min_level = 0):
        # We discard grids if they are ABOVE the level
        mask = na.empty(self.grids.size, dtype='int32')
        get_box_grids_below_level(left_edge, right_edge,
                            level,
                            self.grid_left_edge, self.grid_right_edge,
                            self.grid_levels.astype("int32"), mask, min_level)
        mask = mask.astype("bool")
        return self.grids[mask], na.where(mask)

    def get_periodic_box_grids_below_level(self, left_edge, right_edge, level,
                                           min_level = 0):
        mask = na.zeros(self.grids.shape, dtype='bool')
        dl = self.parameter_file.domain_left_edge
        dr = self.parameter_file.domain_right_edge
        left_edge = na.array(left_edge)
        right_edge = na.array(right_edge)
        dw = dr - dl
        left_dist = left_edge - dl
        db = right_edge - left_edge
        for off_x in [-1, 0, 1]:
            nle = left_edge.copy()
            nle[0] = (dw[0]*off_x + dl[0]) + left_dist[0]
            for off_y in [-1, 0, 1]:
                nle[1] = (dw[1]*off_y + dl[1]) + left_dist[1]
                for off_z in [-1, 0, 1]:
                    nle[2] = (dw[2]*off_z + dl[2]) + left_dist[2]
                    nre = nle + db
                    g, gi = self.get_box_grids_below_level(nle, nre,
                                            level, min_level)
                    mask[gi] = True
        return self.grids[mask], na.where(mask)

