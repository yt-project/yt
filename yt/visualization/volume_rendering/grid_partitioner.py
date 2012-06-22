"""
Import the components of the volume rendering extension

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
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

from yt.utilities.lib import PartitionedGrid, ProtoPrism, GridFace, \
    grid_points_in_volume, find_grids_in_inclined_box
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, parallel_root_only

# HISTORICAL NOTE OF SOME IMPORT:
#   The homogenized brick class should and will be removed.  The more general,
#   and less fancy-parallel, HomogenizedVolume is the way to go.  The HBC is a
#   complicated mechanism that was designed for back when we were going to be
#   passing bricks all around between processors.

class HomogenizedVolume(ParallelAnalysisInterface):
    bricks = None
    def __init__(self, fields = "Density", source = None, pf = None,
                 log_fields = None, no_ghost = False,
                 max_level = 48):
        # Typically, initialized as hanging off a hierarchy.  But, not always.
        ParallelAnalysisInterface.__init__(self)
        self.no_ghost = no_ghost
        if pf is not None: self.pf = pf
        if source is None: source = self.pf.h.all_data()
        self.source = source
        self.fields = ensure_list(fields)
        if log_fields is not None:
            log_fields = ensure_list(log_fields)
        else:
            log_fields = [self.pf.field_info[field].take_log
                         for field in self.fields]
        self.max_level = max_level
        self.log_fields = log_fields

    def traverse(self, back_point, front_point, image):
        mylog.info("Traversing %s bricks between %s and %s",
                   len(self.bricks), back_point, front_point)
        if self.bricks is None: self.initialize_source()
        vec = front_point - back_point
        dist = na.minimum(
             na.sum((self.brick_left_edges - back_point) * vec, axis=1),
             na.sum((self.brick_right_edges - back_point) * vec, axis=1))
        ind = na.argsort(dist)
        for b in self.bricks[ind]:
            #print b.LeftEdge, b.RightEdge
            yield b

    def _partition_grid(self, grid):

        # This is not super efficient, as it re-fills the regions once for each
        # field.
        vcds = []
        for field, log_field in zip(self.fields, self.log_fields):
            vcd = grid.get_vertex_centered_data(field, no_ghost = self.no_ghost)
            vcd = vcd.astype("float64")
            if log_field: vcd = na.log10(vcd)
            vcds.append(vcd)

        GF = GridFaces(grid.Children + [grid])
        PP = ProtoPrism(grid.id, grid.LeftEdge, grid.RightEdge, GF)

        pgs = []
        cm = grid.child_mask.copy()
        if grid.Level > self.max_level:
            return pgs
        elif grid.Level == self.max_level:
            cm[:] = 1
        for P in PP.sweep(0):
            sl = P.get_brick(grid.LeftEdge, grid.dds, cm)
            if len(sl) == 0: continue
            dd = [d[sl[0][0]:sl[0][1]+1,
                    sl[1][0]:sl[1][1]+1,
                    sl[2][0]:sl[2][1]+1].copy() for d in vcds]
            pgs.append(PartitionedGrid(grid.id, len(self.fields), dd,
                        P.LeftEdge, P.RightEdge, sl[-1]))
        return pgs

    def initialize_source(self, source = None):
        if self.bricks is not None and source is not None:
            raise NotImplementedError("Sorry, dynamic shuffling of bricks is" + 
                                      " not yet supported")
        if self.bricks is not None and source is None: return
        bricks = []
        self.comm.preload(self.source._grids, self.fields, self.pf.h.io)
        pbar = get_pbar("Partitioning ", len(self.source._grids))
        for i, g in enumerate(self.source._grids):
            pbar.update(i)
            bricks += self._partition_grid(g)
        pbar.finish()
        self.initialize_bricks(bricks)

    def initialize_bricks(self, bricks):
        NB = len(bricks)
        # Now we set up our (local for now) hierarchy.  Note that to calculate
        # intersection, we only need to do the left edge & right edge.
        #
        # We're going to double up a little bit here in memory.
        self.brick_left_edges = na.zeros( (NB, 3), dtype='float64')
        self.brick_right_edges = na.zeros( (NB, 3), dtype='float64')
        self.brick_parents = na.zeros( NB, dtype='int64')
        self.brick_dimensions = na.zeros( (NB, 3), dtype='int64')
        self.bricks = na.empty(len(bricks), dtype='object')
        for i,b in enumerate(bricks):
            self.brick_left_edges[i,:] = b.LeftEdge
            self.brick_right_edges[i,:] = b.RightEdge
            self.brick_parents[i] = b.parent_grid_id
            self.brick_dimensions[i,:] = b.my_data[0].shape
            self.bricks[i] = b
        # Vertex-centered means we subtract one from the shape
        self.brick_dimensions -= 1

    def reflect_across_boundaries(self):
        mylog.warning("Note that this doesn't fix ghost zones, so there may be artifacts at domain boundaries!")
        nb = []
        # Simplest, clearest iteration ...
        for i in [-1, 1]:
            for j in [-1, 1]:
                for k in [-1, 1]:
                    for b in self.bricks:
                        BB = na.array([b.LeftEdge * [i,j,k], b.RightEdge * [i,j,k]])
                        LE, RE = na.min(BB, axis=0), na.max(BB, axis=0)
                        nb.append(
                            PartitionedGrid(b.parent_grid_id, len(b.my_data), 
                                [md[::i,::j,::k].copy("C") for md in b.my_data],
                                LE, RE, na.array(b.my_data[0].shape) - 1))
        # Replace old bricks
        self.initialize_bricks(nb)

    def store_bricks(self, fn):
        import h5py, cPickle
        f = h5py.File(fn, "w")
        f.create_dataset("/left_edges", data=self.brick_left_edges)
        f.create_dataset("/right_edges", data=self.brick_right_edges)
        f.create_dataset("/parents", data=self.brick_parents)
        f.create_dataset("/dimensions", data=self.brick_dimensions)
        f.create_group("/bricks")
        for i,b in enumerate(self.bricks):
            f.create_group("/bricks/brick_%08i" % i)
            for fi,field in enumerate(self.fields):
                f.create_dataset("/bricks/brick_%08i/%s" % (i, field),
                                 data=b.my_data[fi])
        f.close()

    def load_bricks(self, fn):
        import h5py
        f = h5py.File(fn, "r")
        self.brick_left_edges = f["/left_edges"][:]
        self.brick_right_edges = f["/right_edges"][:]
        self.brick_parents = f["/parents"][:]
        self.brick_dimensions= f["/dimensions"][:]
        bricks = []
        for i,ds in enumerate(sorted(f["/bricks"])):
            td = [f["/bricks/%s/%s" % (ds, field)][:] for field in self.fields]
            bricks.append(PartitionedGrid(
                                self.brick_parents[i], len(td), td,
                                self.brick_left_edges[i,:],
                                self.brick_right_edges[i,:],
                                self.brick_dimensions[i,:],
                                ))
        self.bricks = na.array(bricks, dtype='object')
        f.close()

    def reset_cast(self):
        pass

class SingleBrickVolume(object):
    bricks = None
    def __init__(self, data_array):
        self.bricks = [PartitionedGrid(-1, 1, 
                       [data_array.astype("float64")],
                       na.zeros(3, dtype='float64'),
                       na.ones(3, dtype='float64'),
                       na.array(data_array.shape, dtype='int64')-1)]
        self.brick_dimensions = na.ones((1, 3), dtype='int64')*data_array.shape

    def initialize_source(self):
        pass

    def traverse(self, back, front, image):
        for b in self.bricks: yield b

    def reset_cast(self):
        pass

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

def export_partitioned_grids(grid_list, fn, int_type=na.int64, float_type=na.float64):
    f = h5py.File(fn, "w")
    pbar = get_pbar("Writing Grids", len(grid_list))
    nelem = sum((grid.my_data.size for grid in grid_list))
    ngrids = len(grid_list)
    group = f.create_group("/PGrids")
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
    n_groups = len(f)
    grid_list = []
    dims = f["/PGrids/Dims"][:].astype(int_type)
    left_edges = f["/PGrids/LeftEdges"][:].astype(float_type)
    right_edges = f["/PGrids/RightEdges"][:].astype(float_type)
    data = f["/PGrids/Data"][:].astype(float_type)
    pbar = get_pbar("Reading Grids", dims.shape[0])
    curpos = 0
    for i in xrange(dims.shape[0]):
        gd = dims[i,:]
        gle, gre = left_edges[i,:], right_edges[i,:]
        gdata = data[curpos:curpos+gd.prod()].reshape(gd)
        # Vertex -> Grid, so we -1 from dims in this
        grid_list.append(PartitionedGrid(-1, gdata, gle, gre, gd - 1))
        curpos += gd.prod()
        pbar.update(i)
    pbar.finish()
    f.close()
    return na.array(grid_list, dtype='object')
