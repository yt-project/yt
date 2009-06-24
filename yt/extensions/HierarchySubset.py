"""
A means of extracting a subset of the hierarchy

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Matthew Turk.  All Rights Reserved.

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


from yt.mods import *
import h5py, os.path

import yt.commands as commands

class DummyHierarchy(object):
    pass

class ConstructedRootGrid(object):
    id = -1
    def __init__(self, pf, level, left_edge, right_edge):
        self.pf = pf
        self.hierarchy = DummyHierarchy()
        self.hierarchy.data_style = -1
        self.Level = level
        self.LeftEdge = left_edge
        self.RightEdge = right_edge
        self.index = na.min([grid.get_global_startindex() for grid in
                             pf.h.select_grids(level)], axis=0).astype('int64')
        self.dds = pf.h.select_grids(level)[0].dds.copy()
        dims = (self.RightEdge-self.LeftEdge)/self.dds
        self.ActiveDimensions = dims
        print "Constructing base grid of size %s" % (self.ActiveDimensions)
        self.cg = pf.h.smoothed_covering_grid(level, self.LeftEdge,
                        self.RightEdge, dims=dims)
        self._calculate_child_masks()

    def _calculate_child_masks(self):
        # This might be slow
        grids, grid_ind = self.pf.hierarchy.get_box_grids(
                    self.LeftEdge, self.RightEdge)
        self.Children = [g for g in grids if g.Level == self.Level + 1]
        self.child_mask = na.ones(self.ActiveDimensions, dtype='int32')
        for c in self.Children:
            si = na.maximum(0, na.rint((c.LeftEdge - self.LeftEdge)/self.dds))
            ei = na.minimum(self.ActiveDimensions,
                    na.rint((c.RightEdge - self.LeftEdge)/self.dds))
            self.child_mask[si[0]:ei[0], si[1]:ei[1], si[2]:ei[2]] = 0

    def __getitem__(self, field):
        return self.cg[field]

    def get_global_startindex(self):
        return self.index

    def get_vertex_centered_data(self, field, smoothed=True):
        # We discard the 'smoothed' keyword argument here
        cg = self.pf.h.smoothed_covering_grid(self.Level,
                    self.LeftEdge - self.dds,
                    self.RightEdge + self.dds,
                    dims = self.ActiveDimensions + 2,
                    num_ghost_zones = 1, fields=[field])
        bds = na.array(zip(cg.left_edge+cg.dds/2.0, cg.right_edge-cg.dds/2.0)).ravel()
        interp = lagos.TrilinearFieldInterpolator(cg[field], bds, ['x','y','z'])
        ad = self.ActiveDimensions + 1
        x,y,z = na.mgrid[self.LeftEdge[0]:self.RightEdge[0]:ad[0]*1j,
                         self.LeftEdge[1]:self.RightEdge[1]:ad[1]*1j,
                         self.LeftEdge[2]:self.RightEdge[2]:ad[2]*1j]
        dd = {'x':x,'y':y,'z':z}
        scalars = interp(dd)
        return scalars
        
class ExtractedHierarchy(object):

    def __init__(self, pf, min_level, max_level = -1, offset = None,
                 always_copy=False):
        self.pf = pf
        self.always_copy = always_copy
        self.min_level = min_level
        self.int_offset = na.min([grid.get_global_startindex() for grid in
                             pf.h.select_grids(min_level)], axis=0).astype('float64')
        min_left = na.min([grid.LeftEdge for grid in
                           pf.h.select_grids(min_level)], axis=0).astype('float64')
        max_right = na.max([grid.RightEdge for grid in 
                                   pf.h.select_grids(min_level)], axis=0).astype('float64')
        if offset is None: offset = (max_right + min_left)/2.0
        self.left_edge_offset = offset
        self.mult_factor = 2**min_level
        self.min_left_edge = self._convert_coords(min_left)
        self.max_right_edge = self._convert_coords(max_right)
        if max_level == -1: max_level = pf.h.max_level
        self.max_level = min(max_level, pf.h.max_level)
        self.final_level = self.max_level - self.min_level
        if len(self.pf.h.select_grids(self.min_level)) > 0:
            self._base_grid = ConstructedRootGrid(self.pf, self.min_level,
                               min_left, max_right)
        else: self._base_grid = None
        
    def select_level(self, level):
        if level == 0 and self._base_grid is not None:
            return [self._base_grid]
        return self.pf.h.select_grids(self.min_level + level)

    def get_levels(self):
        for level in range(self.final_level+1):
            yield self.select_level(level)

    def export_output(self, afile, n, field):
        # I prefer dict access, but tables doesn't.
        # But h5py does!
        time_node = afile.create_group("/time-%s" % n)
        time_node.attrs['time'] = self.pf["InitialTime"]
        time_node.attrs['numLevels'] = self.pf.h.max_level+1-self.min_level
        # Can take a while, so let's get a progressbar
        self._export_all_levels(afile, time_node, field)

    def _export_all_levels(self, afile, time_node, field):
        pbar = yt.funcs.get_pbar("Exporting levels", self.final_level+1)
        for i,grid_set in enumerate(self.get_levels()):
            pbar.update(i)
            self.export_level(afile, time_node, i, field, grid_set)
        pbar.finish()

    def export_level(self, afile, time_node, level, field, grids = None):
        level_node = afile.create_group("%s/level-%s" % (time_node,level))
        # Grid objects on this level...
        if grids is None: grids = self.pf.h.select_grids(level+self.min_level)
        level_node.attrs['delta'] = grids[0].dds*self.mult_factor
        level_node.attrs['relativeRefinementFactor'] = na.array([2]*3, dtype='int32')
        level_node.attrs['numGrids'] = len(grids)
        for i,g in enumerate(grids):
            self.export_grid(afile, level_node, g, i, field)

    def _convert_grid(self, grid):
        int_origin = (grid.get_global_startindex() \
                    - self.int_offset*2**(grid.Level-self.min_level)).astype('int64')
        level_int_origin = (grid.LeftEdge - self.left_edge_offset)/grid.dds
        origin = self._convert_coords(grid.LeftEdge)
        dds = grid.dds * self.mult_factor
        return int_origin, level_int_origin, origin, dds

    def export_grid(self, afile, level_node, grid, i, field):
        grid_node = afile.create_group("%s/grid-%s" % (level_node,i))
        int_origin, lint, origin, dds = self._convert_grid(grid)
        grid_node.attrs['integerOrigin'] = int_origin
        grid_node.attrs['origin'] = origin
        grid_node.attrs['ghostzoneFlags'] = na.zeros(6, dtype='int32')
        grid_node.attrs['numGhostzones'] = na.zeros(3, dtype='int32')
        grid_node.attrs['dims'] = grid.ActiveDimensions[::-1].astype('int32')
        if not self.always_copy and self.pf.h.data_style == 6 \
           and field in self.pf.h.field_list:
            if grid.hierarchy.data_style == -1: # constructed grid
                # if we can get conversion in amira we won't need to do this
                ff = grid[field].astype('float32')
                ff /= self.pf.conversion_factors.get(field, 1.0)
                afile.create_dataset("%s/grid-data" % grid_node, data=ff.swapaxes(0,2))
            else:
                tfn = os.path.abspath(afile.filename)
                gfn = os.path.abspath(grid.filename)
                fpn = os.path.commonprefix([tfn, grid.filename])
                fn = grid.filename[len(os.path.commonprefix([tfn, grid.filename])):]
                grid_node.attrs['referenceFileName'] = fn
                grid_node.attrs['referenceDataPath'] = \
                    "/Grid%08i/%s" % (grid.id, field)
        else:
            # Export our array
            afile.create_dataset("%s/grid-data" % grid_node, 
                                 data = grid[field].astype('float32').swapaxes(0,2))

    def _convert_coords(self, val):
        return (val - self.left_edge_offset)*self.mult_factor
