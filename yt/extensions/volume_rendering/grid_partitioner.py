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

from yt.amr_utils import PartitionedGrid, ProtoPrism, GridFace, \
        grid_points_in_volume
from yt.lagos import ParallelAnalysisInterface, only_on_root

class HomogenizedBrickCollection(ParallelAnalysisInterface):
    def __init__(self, source):
        # The idea here is that we have two sources -- the global_domain
        # source, which would be a decomposition of the 3D domain, and a
        # local_domain source, which is the set of bricks we want at the end.
        self.source = source

    @classmethod
    def load_bricks(self, base_filename):
        pass

    def write_my_bricks(self, base_filename):
        pass

    def store_bricks(self, base_filename):
        pass
    
    @parallel_root_only
    def write_hierarchy(self, base_filename):
        pass
    
    def _partition_grid(self, grid, field, log_field = True):
        vcd = grid.get_vertex_centered_data(field).astype('float64')
        if log_field: vcd = na.log10(vcd)

        GF = GridFaces(grid.Children + [grid])
        PP = ProtoPrism(grid.id, grid.LeftEdge, grid.RightEdge, GF)

        pgs = []
        for P in PP.sweep(0):
            pgs += P.get_brick(grid.LeftEdge, grid.dds, vcd, grid.child_mask)
        return pgs

    def _partition_local_grids(self, fields = "Density", log_field = True):
        fields = ensure_list(fields)
        bricks = []
        # We preload.
        self._preload(self._get_grid_objs(), fields, self.pf.h.io)
        pbar = get_pbar("Partitioning ", len(grid_list))
        for i, g in enumerate(self._get_grids()):
            pbar.update(i)
            bricks += self._partition_grid(g, field, log_field, threshold)
        pbar.finish()
        bricks = na.array(bricks, dtype='object')
        NB = len(bricks)
        # Now we set up our (local for now) hierarchy.  Note that to calculate
        # intersection, we only need to do the left edge & right edge.
        self.brick_left_edges = na.fromiter((b.LeftEdge.ravel()
                                             for b in bricks),
                                            dtype='float64', count = NB*3)
        self.brick_left_edges = self.brick_left_edges.reshape((NB,3)
        self.brick_right_edges = na.fromiter((b.RighttEdge.ravel()
                                              for b in bricks),
                                             dtype='float64', count = NB*3)
        self.brick_right_edges = self.brick_right_edges.reshape((NB,3)
        self.brick_parents = na.fromiter((b.parent_grid_id
                                          for b in bricks),
                                         dtype='int64', count = NB)
        self.brick_dimensions = na.fromiter((na.ravel(b.my_data.shape)
                                          for b in bricks),
                                         dtype='int32', count = NB*3)
        self.brick_dimensions = self.brick_dimensions.reshape((NB,3)
        self.brick_owners = na.ones(NB, dtype='int32') * self._mpi_get_rank()
        self.bricks = na.array(bricks, dtype='object')
        self._collect_hierarchy()

    def _get_object_info(self):
        info_dict = dict(left_edges = self.brick_left_edges,
                         right_edges = self.brick_right_edges,
                         parents = self.brick_parents,
                         dims = self.brick_dimensions,)
        return info_dict

    def _set_object_info(self, info_dict):
        self.brick_left_edges = info_dict.pop("left_edges")
        self.brick_right_edges = info_dict.pop("right_edges")
        self.brick_parents = info_dict.pop("parents")
        self._object_parents = self.brick_parents
        self.brick_owners = info_dict.pop("owners")
        bricks = self.bricks
        self.bricks = na.array([None] * self.brick_owners.size, dtype='object')
        # Copy our bricks back in
        self.bricks[self.brick_owners = self._mpi_get_rank()] = bricks[:]

    def _pack_object(self, index):
        return self.bricks[index].my_data.ravel()

    def _unpack_object(self, index, data):
        pgi = self.brick_parents[index]
        LE = self.brick_left_edges[index,:].copy()
        RE = self.brick_right_edges[index,:].copy()
        dims = self.brick_dimensions[index,:].copy()
        data = data.reshape(dims)
        self.bricks[index] = PartitionedGrid(
                pgi, data, LE, RE, dims)

    def _wipe_objects(self, indices):
        self.bricks[indices] = None

    def _collect_bricks(self, intersection_source):
        if not self._distributed: return
        # This entire routine should instead be set up to do:
        #   alltoall broadcast of the *number* of requested bricks
        #   non-blocking receives posted for int arrays
        #   sizes of data calculated
        #   concatenated data receives posted
        #   send all data
        #   get bricks back
        # This presupposes that we are using the AMRInclinedBox as a data
        # source.  If we're not, we ought to be.
        needed_brick_i = find_grids_in_inclined_box(
            intersection_source.box_vectors, intersection_source.center,
            self.brick_left_edges, self.brick_right_edges)
        needed_brick_i = na.where(needed_brick_i)[0]
        self._collect_objects(needed_brick_i)

    def _initialize_parallel(self):
        pass

    def _finalize_parallel(self):
        pass

    def get_brick(self, brick_id):
        pass

    @property
    def _grids(self):
        return self.source._grids

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
    n_groups = len(f.listnames())
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
        rid_list.append(PartitionedGrid(gdata, gle, gre, gd - 1))
        curpos += gd.prod()
        pbar.update(i)
    pbar.finish()
    f.close()
    return na.array(grid_list, dtype='object')
