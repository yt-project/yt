"""
Import the components of the volume rendering extension

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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
import h5py

from yt.amr_utils import PartitionedGrid, ProtoPrism, GridFace

class GridFaces(object):
    def __init__(self, grids):
        self.faces = [ [], [], [] ]
        for grid in grids:
            for direction in range(3):
                self.faces[direction].append( GridFace(grid, direction, 1) )
                self.faces[direction].append( GridFace(grid, direction, 0) )
        for f in self.faces:
            f.sort(key = lambda a: a.coord)

    def __getitem__(self, item):
        return self.faces[item]

def partition_grid(start_grid, field, log_field = True, threshold = None):
    if threshold is not None:
        if start_grid[field].max() < threshold[0] or \
           start_grid[field].min() > threshold[1]: return None
    to_cut_up = start_grid.get_vertex_centered_data(field, smoothed=True).astype('float64')

    if log_field: to_cut_up = na.log10(to_cut_up)

    GF = GridFaces(start_grid.Children + [start_grid])
    PP = ProtoPrism(start_grid.LeftEdge, start_grid.RightEdge, GF)
    pgs = []
    for P in PP.sweep(0):
        pgs += P.get_brick(start_grid.LeftEdge, start_grid.dds, to_cut_up, start_grid.child_mask)
    return pgs

    if len(start_grid.Children) == 0:
        pg = PartitionedGrid(
                to_cut_up.copy(),
                na.array(start_grid.LeftEdge, dtype='float64'),
                na.array(start_grid.RightEdge, dtype='float64'),
                na.array(start_grid.ActiveDimensions, dtype='int64'))
        return [pg]

    x_vert = [0, start_grid.ActiveDimensions[0]]
    y_vert = [0, start_grid.ActiveDimensions[1]]
    z_vert = [0, start_grid.ActiveDimensions[2]]

    gi = start_grid.get_global_startindex()
    for grid in start_grid.Children:
        si = grid.get_global_startindex()/2 - gi
        ei = si + grid.ActiveDimensions/2 
        x_vert += [si[0], ei[0]]
        y_vert += [si[1], ei[1]]
        z_vert += [si[2], ei[2]]

    # Now we sort by our vertices, in axis order

    x_vert.sort()
    y_vert.sort()
    z_vert.sort()

    return [g for g in _partition(start_grid, to_cut_up, x_vert, y_vert, z_vert)]

def _partition(grid, grid_data, x_vert, y_vert, z_vert):
    grids = []
    cim = grid.child_index_mask
    for xs, xe in zip(x_vert[:-1], x_vert[1:]):
        for ys, ye in zip(y_vert[:-1], y_vert[1:]):
            for zs, ze in zip(z_vert[:-1], z_vert[1:]):
                sl = (slice(xs, xe), slice(ys, ye), slice(zs, ze))
                dd = cim[sl]
                if dd.size == 0: continue
                uniq = na.unique(dd)
                if uniq.size > 1: continue
                if uniq[0] > -1: continue
                data = grid_data[xs:xe+1,ys:ye+1,zs:ze+1].copy()
                dims = na.array(dd.shape, dtype='int64')
                start_index = na.array([xs,ys,zs], dtype='int64')
                left_edge = grid.LeftEdge + start_index * grid.dds
                right_edge = left_edge + dims * grid.dds
                yield PartitionedGrid(
                    data, left_edge, right_edge, dims)

def partition_all_grids(grid_list, field = "Density", log_field = True,
                        threshold = (-1e300, 1e300), eval_func = None):
    new_grids = []
    pbar = get_pbar("Partitioning ", len(grid_list))
    if eval_func is None: eval_func = lambda a: True
    dx = 1e300
    for i, g in enumerate(grid_list):
        if not eval_func(g): continue
        pbar.update(i)
        if g.dds[0] < dx: dx = g.dds[0]
        to_add = partition_grid(g, field, log_field, threshold)
        if to_add is not None: new_grids += to_add
    pbar.finish()
    for g in new_grids: g.min_dds = dx
    return na.array(new_grids, dtype='object')

def export_partitioned_grids(grid_list, fn, int_type=na.int64, float_type=na.float64):
    f = h5py.File(fn, "w")
    pbar = get_pbar("Writing Grids", len(grid_list))
    nelem = sum((grid.my_data.size for grid in grid_list))
    ngrids = len(grid_list)
    group = f.create_group("/PGrids")
    group.attrs["min_dds"] = grid_list[0].min_dds
    left_edge = na.concatenate([[grid.LeftEdge,] for grid in grid_list])
    f.create_dataset("/PGrids/LeftEdges", data=left_edge, dtype=float_type); del left_edge
    right_edge = na.concatenate([[grid.RightEdge,] for grid in grid_list])
    f.create_dataset("/PGrids/RightEdges", data=right_edge, dtype=float_type); del right_edge
    dims = na.concatenate([[grid.my_data.shape[:],] for grid in grid_list])
    f.create_dataset("/PGrids/Dims", data=dims, dtype=int_type); del dims
    data = na.concatenate([grid.my_data.ravel() for grid in grid_list])
    f.create_dataset("/PGrids/Data", data=data, dtype=float_type); del data
    f.close()
    pbar.finish()

def import_partitioned_grids(fn, int_type=na.int64, float_type=na.float64):
    f = h5py.File(fn, "r")
    n_groups = len(f.listnames())
    grid_list = []
    dims = f["/PGrids/Dims"][:].astype(int_type)
    left_edges = f["/PGrids/LeftEdges"][:].astype(float_type)
    right_edges = f["/PGrids/RightEdges"][:].astype(float_type)
    data = f["/PGrids/Data"][:].astype(float_type)
    pbar = get_pbar("Reading Grids", dims.shape[0])
    curpos = 0
    dx = f["/PGrids"].attrs["min_dds"]
    for i in xrange(dims.shape[0]):
        gd = dims[i,:]
        gle, gre = left_edges[i,:], right_edges[i,:]
        gdata = data[curpos:curpos+gd.prod()].reshape(gd)
        # Vertex -> Grid, so we -1 from dims in this
        grid_list.append(PartitionedGrid(gdata, gle, gre, gd - 1))
        grid_list[-1].min_dds = dx
        curpos += gd.prod()
        pbar.update(i)
    pbar.finish()
    f.close()
    return na.array(grid_list, dtype='object')

class PartitionRegion(object):
    _count = 0
    def __init__(self, dims, source_offset, source_vertices, cim_base):
        self.source_offset = source_offset
        self.dims = dims
        cv = []
        self._cim = cim_base
        self.child_vertices = source_vertices

    @property
    def cim(self):
        return self._cim[self.sl]

    @property
    def sl(self):
        sls = self.source_offset
        sle = self.source_offset + self.dims
        return tuple([slice(sls[i], sle[i]) for i in range(3)])

    def split(self, axis, coord):
        dims_left = self.dims.copy()
        dims_left[axis] = coord - self.source_offset[axis]
        off_left = self.source_offset.copy()
        left_region = PartitionRegion(dims_left, off_left,
                        self.child_vertices, self._cim)
        dims_right = self.dims.copy()
        dims_right[axis] = self.dims[axis] - coord + self.source_offset[axis]
        off_right = self.source_offset.copy()
        off_right[axis] = coord
        right_region = PartitionRegion(dims_right, off_right,
                        self.child_vertices, self._cim)
        return left_region, right_region
        
    def find_hyperplane(self, axis):
        # Our axis is the normal to the hyperplane
        # Region boundaries is [2][3]
        # child_vertices is flat 3D array
        min_balance = 1e30
        considered = set([self.source_offset[axis]])
        considered.add(self.source_offset[axis] + self.dims[axis])
        best_coord = self.source_offset[axis] + self.dims[axis]
        for v in self.child_vertices:
            coord = v[axis]
            sc = coord - self.source_offset[axis]
            if coord in considered: continue
            if sc >= self.dims[axis]: continue
            if sc < 0: continue
            eff = self.evaluate_hyperplane(axis, coord)
            if eff < min_balance:
                min_balance = eff
                best_coord = coord
            considered.add(coord)
        return best_coord

    def evaluate_hyperplane(self, axis, coord):
        # We check that we're roughly evenly balanced on either side of the grid
        # Calculate how well balanced it is...
        vert = self.child_vertices[:,axis]
        n_left = (vert <= coord).sum()
        n_right = (vert > coord).sum()
        eff = abs(0.5 - (n_left / float(vert.shape[0])))
        return eff

def partition_region(region, axis=0):
    # region_boundaries is in ints
    split_coord = region.find_hyperplane(axis)
    sc = split_coord - region.source_offset[axis]
    if sc == 0 or sc == region.dims[axis]:
        rc = na.unique(region.cim)
        if rc.size > 1 and rc[0] == -1:
            region._count += 1
            if region._count > 3:
                import pdb;pdb.set_trace()
            return partition_region(region, (axis+1)%3)
        elif rc.size > 1 and rc[0] > -1:
            return []
    left_region, right_region = region.split(axis, split_coord)
    lrc = na.unique(left_region.cim)
    rrc = na.unique(right_region.cim)
    if lrc.size > 1:
        if lrc[0] == -1:
            left_region = partition_region(left_region, (axis + 1) % 3)
        if lrc[0] > -1:
            left_region = []
        #print axis, split_coord, "Not splitting left region", lrc
    else:
        if lrc[0] == -1:
            left_region = [left_region]
        else:
            left_region = []

    if rrc.size > 1:
        if rrc[0] == -1:
            right_region = partition_region(right_region, (axis + 1) % 3)
        if rrc[0] > -1:
            right_region = []
    else:
        if rrc[0] == -1:
            right_region = [right_region]
        else:
            right_region = []

    return left_region + right_region
