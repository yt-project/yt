"""
AMR index container class



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.config import ytcfg
from yt.funcs import \
    mylog, ensure_numpy_array, \
    ensure_list
from yt.utilities.lib.misc_utilities import \
    get_box_grids_below_level
from yt.geometry.grid_container import \
    MatchPointsToGrids, \
    GridTree
from yt.utilities.physical_ratios import \
    HUGE
from yt.utilities.exceptions import YTTooParallel

class ObjectFindingMixin(object) :

    def find_ray_grids(self, coord, axis):
        """
        Returns the (objects, indices) of grids that an (x,y) ray intersects
        along *axis*
        """
        # Let's figure out which grids are on the slice
        mask=np.ones(self.num_grids)
        # So if gRE > coord, we get a mask, if not, we get a zero
        #    if gLE > coord, we get a zero, if not, mask
        # Thus, if the coordinate is between the two edges, we win!
        xax = self.ds.coordinates.x_axis[axis]
        yax = self.ds.coordinates.y_axis[axis]
        np.choose(np.greater(self.grid_right_edge[:,xax],coord[0]),(0,mask),mask)
        np.choose(np.greater(self.grid_left_edge[:,xax],coord[0]),(mask,0),mask)
        np.choose(np.greater(self.grid_right_edge[:,yax],coord[1]),(0,mask),mask)
        np.choose(np.greater(self.grid_left_edge[:,yax],coord[1]),(mask,0),mask)
        ind = np.where(mask == 1)
        return self.grids[ind], ind

    def find_max(self, field, finest_levels = 3):
        """
        Returns (value, center) of location of maximum for a given field.
        """
        if (field, finest_levels) in self._max_locations:
            return self._max_locations[(field, finest_levels)]
        mv, pos = self.find_max_cell_location(field, finest_levels)
        self._max_locations[(field, finest_levels)] = (mv, pos)
        return mv, pos

    def find_max_cell_location(self, field, finest_levels = 3):
        if finest_levels is not False:
            # This prevents bad values for the case that the number of grids to
            # search is smaller than the number of processors being applied to
            # the task, by 
            nproc = ytcfg.getint("yt", "__topcomm_parallel_size")
            while 1:
                gi = (self.grid_levels >= self.max_level - finest_levels).ravel()
                if gi.sum() >= nproc:
                    break
                elif finest_levels >= self.max_level:
                    raise YTTooParallel
                else:
                    finest_levels += 1
                
            source = self.grid_collection([0.0]*3, self.grids[gi])
        else:
            source = self.all_data()
        mylog.debug("Searching %s grids for maximum value of %s",
                    len(source._grids), field)
        max_val, mx, my, mz = \
            source.quantities["MaxLocation"]( field )
        mylog.info("Max Value is %0.5e at %0.16f %0.16f %0.16f", 
              max_val, mx, my, mz)
        self.parameters["Max%sValue" % (field)] = max_val
        self.parameters["Max%sPos" % (field)] = "%s" % ((mx,my,mz),)
        return max_val, np.array((mx,my,mz), dtype='float64')

    def find_min(self, field):
        """
        Returns (value, center) of location of minimum for a given field
        """
        gI = np.where(self.grid_levels >= 0) # Slow but pedantic
        minVal = HUGE
        for grid in self.grids[gI[0]]:
            mylog.debug("Checking %s (level %s)", grid.id, grid.Level)
            val, coord = grid.find_min(field)
            if val < minVal:
                minCoord = coord
                minVal = val
                minGrid = grid
        mc = np.array(minCoord)
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
        mask=np.ones(self.num_grids)
        for i in range(len(coord)):
            np.choose(np.greater(self.grid_left_edge[:,i],coord[i]), (mask,0), mask)
            np.choose(np.greater(self.grid_right_edge[:,i],coord[i]), (0,mask), mask)
        ind = np.where(mask == 1)
        return self.grids[ind], ind

    def find_points(self, x, y, z) :
        """
        Returns the (objects, indices) of leaf grids containing a number of (x,y,z) points
        """
        x = ensure_numpy_array(x)
        y = ensure_numpy_array(y)
        z = ensure_numpy_array(z)
        if not len(x) == len(y) == len(z):
            raise AssertionError("Arrays of indices must be of the same size")

        grid_tree = self.get_grid_tree()
        pts = MatchPointsToGrids(grid_tree, len(x), x, y, z)
        ind = pts.find_points_in_tree()
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
        >>> ds.h.find_field_value_at_point(['Density', 'Temperature'],
            [0.4, 0.3, 0.8])
        [2.1489e-24, 1.23843e4]
        """
        # Get the most-refined grid at this coordinate.
        this = self.find_point(coord)[0][-1]
        cellwidth = (this.RightEdge - this.LeftEdge) / this.ActiveDimensions
        mark = np.zeros(3).astype('int')
        # Find the index for the cell containing this point.
        for dim in range(len(coord)):
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
        mask = np.ones(self.num_grids)
        # So if gRE > coord, we get a mask, if not, we get a zero
        #    if gLE > coord, we get a zero, if not, mask
        # Thus, if the coordinate is between the edges, we win!
        np.choose(np.greater(self.grid_right_edge[:,axis],coord),(0,mask),mask)
        np.choose(np.greater(self.grid_left_edge[:,axis],coord),(mask,0),mask)
        ind = np.where(mask == 1)
        return self.grids[ind], ind

    def find_sphere_grids(self, center, radius):
        """
        Returns objects, indices of grids within a sphere
        """
        centers = (self.grid_right_edge + self.grid_left_edge)/2.0
        long_axis = np.maximum.reduce(self.grid_right_edge - self.grid_left_edge, 1)
        t = np.abs(centers - center)
        DW = self.dataset.domain_right_edge \
           - self.dataset.domain_left_edge
        np.minimum(t, np.abs(DW-t), t)
        dist = np.sqrt(np.sum((t**2.0), axis=1))
        gridI = np.where(dist < (radius + long_axis))
        return self.grids[gridI], gridI

    def get_box_grids(self, left_edge, right_edge):
        """
        Gets back all the grids between a left edge and right edge
        """
        eps = np.finfo(np.float64).eps
        grid_i = np.where(
            (np.all((self.grid_right_edge - left_edge) > eps, axis=1) &
             np.all((right_edge - self.grid_left_edge) > eps, axis=1))
        )

        return self.grids[grid_i], grid_i

    def get_periodic_box_grids(self, left_edge, right_edge):
        mask = np.zeros(self.grids.shape, dtype='bool')
        dl = self.dataset.domain_left_edge
        dr = self.dataset.domain_right_edge
        left_edge = np.array(left_edge)
        right_edge = np.array(right_edge)
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
        return self.grids[mask], np.where(mask)

    def get_box_grids_below_level(self, left_edge, right_edge, level,
                                  min_level = 0):
        # We discard grids if they are ABOVE the level
        mask = np.empty(self.grids.size, dtype='int32')
        get_box_grids_below_level(left_edge, right_edge,
                            level,
                            self.grid_left_edge, self.grid_right_edge,
                            self.grid_levels.astype("int32"), mask, min_level)
        mask = mask.astype("bool")
        return self.grids[mask], np.where(mask)

    def get_periodic_box_grids_below_level(self, left_edge, right_edge, level,
                                           min_level = 0):
        mask = np.zeros(self.grids.shape, dtype='bool')
        dl = self.dataset.domain_left_edge
        dr = self.dataset.domain_right_edge
        left_edge = np.array(left_edge)
        right_edge = np.array(right_edge)
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
        return self.grids[mask], np.where(mask)

    def get_grid_tree(self) :

        left_edge = np.zeros((self.num_grids, 3))
        right_edge = np.zeros((self.num_grids, 3))
        level = np.zeros((self.num_grids), dtype='int64')
        parent_ind = np.zeros((self.num_grids), dtype='int64')
        num_children = np.zeros((self.num_grids), dtype='int64')
        dimensions = np.zeros((self.num_grids, 3), dtype="int32")

        for i, grid in enumerate(self.grids) :

            left_edge[i,:] = grid.LeftEdge
            right_edge[i,:] = grid.RightEdge
            level[i] = grid.Level
            if grid.Parent is None :
                parent_ind[i] = -1
            else :
                parent_ind[i] = grid.Parent.id - grid.Parent._id_offset
            num_children[i] = np.int64(len(grid.Children))
            dimensions[i,:] = grid.ActiveDimensions

        return GridTree(self.num_grids, left_edge, right_edge, dimensions,
                        parent_ind, level, num_children)
