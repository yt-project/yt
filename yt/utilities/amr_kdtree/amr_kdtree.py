"""
AMR kD-Tree Framework

Authors: Samuel Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder

Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Samuel Skillman.  All Rights Reserved.

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
from yt.funcs import *
import numpy as np
import h5py
from amr_kdtools import Node, Split, kd_is_leaf, kd_sum_volume, kd_node_check, \
        depth_traverse, viewpoint_traverse, add_grids, \
        receive_and_reduce, send_to_parent, scatter_image, find_node, \
        depth_first_touch
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import ParallelAnalysisInterface 
from yt.utilities.lib.grid_traversal import PartitionedGrid
from yt.utilities.math_utils import periodic_position

steps = np.array([[-1, -1, -1], [-1, -1,  0], [-1, -1,  1],
                  [-1,  0, -1], [-1,  0,  0], [-1,  0,  1],
                  [-1,  1, -1], [-1,  1,  0], [-1,  1,  1],
                  
                  [ 0, -1, -1], [ 0, -1,  0], [ 0, -1,  1],
                  [ 0,  0, -1],
                  # [ 0,  0,  0],
                  [ 0,  0,  1],
                  [ 0,  1, -1], [ 0,  1,  0], [ 0,  1,  1],
                  
                  [ 1, -1, -1], [ 1, -1,  0], [ 1, -1,  1],
                  [ 1,  0, -1], [ 1,  0,  0], [ 1,  0,  1],
                  [ 1,  1, -1], [ 1,  1,  0], [ 1,  1,  1] ])


def make_vcd(data):
    new_field = np.zeros(np.array(data.shape) + 1, dtype='float64')
    of = data
    new_field[:-1, :-1, :-1] += of
    new_field[:-1, :-1, 1:] += of
    new_field[:-1, 1:, :-1] += of
    new_field[:-1, 1:, 1:] += of
    new_field[1:, :-1, :-1] += of
    new_field[1:, :-1, 1:] += of
    new_field[1:, 1:, :-1] += of
    new_field[1:, 1:, 1:] += of
    np.multiply(new_field, 0.125, new_field)

    new_field[:, :, -1] = 2.0*new_field[:, :, -2] - new_field[:, :, -3]
    new_field[:, :, 0] = 2.0*new_field[:, :, 1] - new_field[:, :, 2]
    new_field[:, -1, :] = 2.0*new_field[:, -2, :] - new_field[:, -3, :]
    new_field[:, 0, :] = 2.0*new_field[:, 1, :] - new_field[:, 2, :]
    new_field[-1, :, :] = 2.0*new_field[-2, :, :] - new_field[-3, :, :]
    new_field[0, :, :] = 2.0*new_field[1, :, :] - new_field[2, :, :]
    return new_field

class Tree(object):
    def __init__(self, pf, comm_rank=0, comm_size=1,
            min_level=None, max_level=None, source=None):
        
        self.pf = pf
        if source is None:
            source = pf.h.all_data()
        self.source = source
        self._id_offset = self.pf.h.grids[0]._id_offset
        if min_level is None: min_level = 0
        if max_level is None: max_level = pf.h.max_level
        self.min_level = min_level
        self.max_level = max_level
        self.comm_rank = comm_rank
        self.comm_size = comm_size
        left_edge = self.source.left_edge
        right_edge= self.source.right_edge
        self.trunk = Node(None, None, None,
                left_edge, right_edge, None, 1)
        self.build()

    def add_grids(self, grids):
        gles = np.array([g.LeftEdge for g in grids])
        gres = np.array([g.RightEdge for g in grids])
        gids = np.array([g.id for g in grids])
        add_grids(self.trunk, gles, gres, gids, self.comm_rank, self.comm_size)
        del gles, gres, gids, grids

    def build(self):
        lvl_range = range(self.min_level, self.max_level+1)
        for lvl in lvl_range:
            #grids = self.source.select_grids(lvl)
            grids = np.array([b for b, mask in self.source.blocks])
            if len(grids) == 0: break
            self.add_grids(grids)

    def check_tree(self):
        for node in depth_traverse(self):
            if node.grid is None:
                continue
            grid = self.pf.h.grids[node.grid - self._id_offset]
            dds = grid.dds
            gle = grid.LeftEdge
            gre = grid.RightEdge
            li = np.rint((node.left_edge-gle)/dds).astype('int32')
            ri = np.rint((node.right_edge-gle)/dds).astype('int32')
            dims = (ri - li).astype('int32')
            assert(np.all(grid.LeftEdge <= node.left_edge))
            assert(np.all(grid.RightEdge >= node.right_edge))
            assert(np.all(dims > 0))
            # print grid, dims, li, ri

        # Calculate the Volume
        vol = kd_sum_volume(self.trunk)
        mylog.debug('AMRKDTree volume = %e' % vol)
        kd_node_check(self.trunk)

    def sum_cells(self):
        cells = 0
        for node in depth_traverse(self):
            if node.grid is None:
                continue
            grid = self.pf.h.grids[node.grid - self._id_offset]
            dds = grid.dds
            gle = grid.LeftEdge
            gre = grid.RightEdge
            li = np.rint((node.left_edge-gle)/dds).astype('int32')
            ri = np.rint((node.right_edge-gle)/dds).astype('int32')
            dims = (ri - li).astype('int32')
            cells += np.prod(dims)

        return cells

class AMRKDTree(ParallelAnalysisInterface):
    fields = None
    log_fields = None
    no_ghost = True
    def __init__(self, pf, min_level=None, max_level=None, source=None):

        ParallelAnalysisInterface.__init__(self)

        self.pf = pf
        self.current_vcds = []
        self.current_saved_grids = []
        self.bricks = []
        self.brick_dimensions = []
        self.sdx = pf.h.get_smallest_dx()
        self._initialized = False
        self._id_offset = pf.h.grids[0]._id_offset

        #self.add_mask_field()
        if source is None:
            source = pf.h.all_data()
        self.source = source
    
        mylog.debug('Building AMRKDTree')
        self.tree = Tree(pf, self.comm.rank, self.comm.size,
                         min_level=min_level,
                         max_level=max_level, source=source)

    def set_fields(self, fields, log_fields, no_ghost):
        self.fields = fields
        self.log_fields = log_fields
        self.no_ghost = no_ghost
        del self.bricks, self.brick_dimensions
        self.brick_dimensions = []
        bricks = []
        for b in self.traverse():
            bricks.append(b)
        self.bricks = np.array(bricks)
        self.brick_dimensions = np.array(self.brick_dimensions)
    
    def initialize_source(self, fields, log_fields, no_ghost):
        if fields == self.fields and log_fields == self.log_fields and \
            no_ghost == self.no_ghost:
            return
        self.set_fields(fields, log_fields, no_ghost)

    def traverse(self, viewpoint=None):
        if viewpoint is None:
            for node in depth_traverse(self.tree):
                if kd_is_leaf(node) and node.grid is not None:
                    yield self.get_brick_data(node)
        else:
            for node in viewpoint_traverse(self.tree, viewpoint):
                if kd_is_leaf(node) and node.grid is not None:
                    yield self.get_brick_data(node)

    def get_node(self, nodeid):
        path = np.binary_repr(nodeid)
        depth = 1
        temp = self.tree.trunk
        for depth in range(1,len(path)):
            if path[depth] == '0':
                temp = temp.left
            else:
                temp = temp.right
        assert(temp is not None)
        return temp

    def locate_node(self, pos):
        return find_node(self.tree.trunk, pos)

    def get_reduce_owners(self):
        owners = {}
        for bottom_id in range(self.comm.size, 2*self.comm.size):
            temp = self.get_node(bottom_id)
            owners[temp.id] = temp.id - self.comm.size
            while temp is not None:
                if temp.parent is None: break
                if temp == temp.parent.right:
                    break
                temp = temp.parent
                owners[temp.id] = owners[temp.left.id]
        return owners

    def reduce_tree_images(self, image, viewpoint):
        if self.comm.size <= 1: return image
        myrank = self.comm.rank
        nprocs = self.comm.size
        owners = self.get_reduce_owners()
        node = self.get_node(nprocs + myrank)

        while True:
            if owners[node.parent.id] == myrank:
                split = node.parent.split
                left_in_front = viewpoint[split.dim] < node.parent.split.pos
                #add_to_front = (left_in_front == (node == node.parent.right))
                add_to_front = not left_in_front
                image = receive_and_reduce(self.comm, owners[node.parent.right.id],
                                  image, add_to_front)
                if node.parent.id == 1: break
                else: node = node.parent
            else:
                send_to_parent(self.comm, owners[node.parent.id], image)
                break
        image = scatter_image(self.comm, owners[1], image)
        return image

    def get_brick_data(self, node):
        if node.data is not None: return node.data
        grid = self.pf.h.grids[node.grid - self._id_offset]
        dds = grid.dds
        gle = grid.LeftEdge
        gre = grid.RightEdge
        li = np.rint((node.left_edge-gle)/dds).astype('int32')
        ri = np.rint((node.right_edge-gle)/dds).astype('int32')
        dims = (ri - li).astype('int32')
        assert(np.all(grid.LeftEdge <= node.left_edge))
        assert(np.all(grid.RightEdge >= node.right_edge))

        if grid in self.current_saved_grids:
            dds = self.current_vcds[self.current_saved_grids.index(grid)]
        else:
            dds = []
            mask = make_vcd(grid.child_mask)
            mask = np.clip(mask, 0.0, 1.0)
            mask[mask<1.0] = np.inf
            for i,field in enumerate(self.fields):
                vcd = make_vcd(grid[field])
                vcd *= mask
                if self.log_fields[i]: vcd = np.log10(vcd)
                dds.append(vcd)
                self.current_saved_grids.append(grid)
                self.current_vcds.append(dds)

        data = [d[li[0]:ri[0]+1,
                  li[1]:ri[1]+1,
                  li[2]:ri[2]+1].copy() for d in dds]

        brick = PartitionedGrid(grid.id, data,
                                node.left_edge.copy(),
                                node.right_edge.copy(),
                                dims.astype('int64'))
        node.data = brick
        if not self._initialized: 
            self.brick_dimensions.append(dims)
        return brick

    def locate_brick(self, position):
        r"""Given a position, find the node that contains it.
        Alias of AMRKDTree.locate_node, to preserve backwards
        compatibility.
        """
        return self.locate_node(position) 

    def locate_neighbors(self, grid, ci):
        r"""Given a grid and cell index, finds the 26 neighbor grids 
        and cell indices.
        
        Parameters
        ----------
        grid: Grid Object
            Grid containing the cell of interest
        ci: array-like
            The cell index of the cell of interest

        Returns
        -------
        grids: Numpy array of Grid objects
        cis: List of neighbor cell index tuples

        Both of these are neighbors that, relative to the current cell
        index (i,j,k), are ordered as: 
        
        (i-1, j-1, k-1), (i-1, j-1, k ), (i-1, j-1, k+1), ...  
        (i-1, j  , k-1), (i-1, j  , k ), (i-1, j  , k+1), ...  
        (i+1, j+1, k-1), (i-1, j-1, k ), (i+1, j+1, k+1)

        That is they start from the lower left and proceed to upper
        right varying the third index most frequently. Note that the
        center cell (i,j,k) is ommitted.
        
        """
        ci = np.array(ci)
        center_dds = grid.dds
        position = grid.LeftEdge + (np.array(ci)+0.5)*grid.dds
        grids = np.empty(26, dtype='object')
        cis = np.empty([26,3], dtype='int64')
        offs = 0.5*(center_dds + self.sdx)

        new_cis = ci + steps
        in_grid = np.all((new_cis >=0)*
                         (new_cis < grid.ActiveDimensions),axis=1)
        new_positions = position + steps*offs
        new_positions = [periodic_position(p, self.pf) for p in new_positions]
        grids[in_grid] = grid
                
        get_them = np.argwhere(in_grid != True).ravel()
        cis[in_grid] = new_cis[in_grid]

        if (in_grid != True).sum()>0:
            grids[in_grid != True] = \
                [self.pf.h.grids[self.locate_brick(new_positions[i]).grid] 
                 for i in get_them]
            cis[in_grid != True] = \
                [(new_positions[i]-grids[i].LeftEdge)/
                 grids[i].dds for i in get_them]
        cis = [tuple(ci) for ci in cis]
        return grids, cis

    def locate_neighbors_from_position(self, position):
        r"""Given a position, finds the 26 neighbor grids 
        and cell indices.

        This is a mostly a wrapper for locate_neighbors.
        
        Parameters
        ----------
        position: array-like
            Position of interest

        Returns
        -------
        grids: Numpy array of Grid objects
        cis: List of neighbor cell index tuples

        Both of these are neighbors that, relative to the current cell
        index (i,j,k), are ordered as: 
        
        (i-1, j-1, k-1), (i-1, j-1, k ), (i-1, j-1, k+1), ...  
        (i-1, j  , k-1), (i-1, j  , k ), (i-1, j  , k+1), ...  
        (i+1, j+1, k-1), (i-1, j-1, k ), (i+1, j+1, k+1)

        That is they start from the lower left and proceed to upper
        right varying the third index most frequently. Note that the
        center cell (i,j,k) is ommitted.
        
        """
        position = np.array(position)
        grid = self.pf.h.grids[self.locate_brick(position).grid]
        ci = ((position-grid.LeftEdge)/grid.dds).astype('int64')
        return self.locate_neighbors(grid,ci)

    def store_kd_bricks(self, fn=None):
        if not self._initialized:
            self.initialize_source()
        if fn is None:
            fn = '%s_kd_bricks.h5'%self.pf
        if self.comm.rank != 0:
            self.comm.recv_array(self.comm.rank-1, tag=self.comm.rank-1)
        f = h5py.File(fn,'w')
        for node in depth_traverse(self.tree):
            i = node.id
            if node.data is not None:
                for fi,field in enumerate(self.fields):
                    try:
                        f.create_dataset("/brick_%s_%s" % (hex(i),field),
                                         data = node.data.my_data[fi].astype('float64'))
                    except:
                        pass
        f.close()
        del f
        if self.comm.rank != (self.comm.size-1):
            self.comm.send_array([0],self.comm.rank+1, tag=self.comm.rank)
        
    def load_kd_bricks(self,fn=None):
        if fn is None:
            fn = '%s_kd_bricks.h5' % self.pf
        if self.comm.rank != 0:
            self.comm.recv_array(self.comm.rank-1, tag=self.comm.rank-1)
        try:
            f = h5py.File(fn,"a")
            for node in depth_traverse(self.tree):
                i = node.id
                if node.grid is not None:
                    data = [f["brick_%s_%s" %
                              (hex(i), field)][:].astype('float64') for field in self.fields]
                    node.data = PartitionedGrid(node.grid.id, data,
                                                 node.l_corner.copy(), 
                                                 node.r_corner.copy(), 
                                                 node.dims.astype('int64'))
                    
                    self.bricks.append(node.data)
                    self.brick_dimensions.append(node.dims)

            self.bricks = np.array(self.bricks)
            self.brick_dimensions = np.array(self.brick_dimensions)

            self._initialized=True
            f.close()
            del f
        except:
            pass
        if self.comm.rank != (self.comm.size-1):
            self.comm.send_array([0],self.comm.rank+1, tag=self.comm.rank)

    def join_parallel_trees(self):
        if self.comm.size == 0: return
        nid, pid, lid, rid, les, res, gid, splitdims, splitposs = \
                self.get_node_arrays()
        nid = self.comm.par_combine_object(nid, 'cat', 'list') 
        pid = self.comm.par_combine_object(pid, 'cat', 'list') 
        lid = self.comm.par_combine_object(lid, 'cat', 'list') 
        rid = self.comm.par_combine_object(rid, 'cat', 'list') 
        gid = self.comm.par_combine_object(gid, 'cat', 'list') 
        les = self.comm.par_combine_object(les, 'cat', 'list') 
        res = self.comm.par_combine_object(res, 'cat', 'list') 
        splitdims = self.comm.par_combine_object(splitdims, 'cat', 'list') 
        splitposs = self.comm.par_combine_object(splitposs, 'cat', 'list') 
        nid = np.array(nid)
        new_tree = self.rebuild_tree_from_array(nid, pid, lid, 
            rid, les, res, gid, splitdims, splitposs)

    def get_node_arrays(self):
        nids = []
        leftids = []
        rightids = []
        parentids = []
        les = []
        res = []
        gridids = []
        splitdims = []
        splitposs = []
        for node in depth_first_touch(self.tree):
            nids.append(node.id) 
            les.append(node.left_edge) 
            res.append(node.right_edge) 
            if node.left is None:
                leftids.append(-1) 
            else:
                leftids.append(node.left.id) 
            if node.right is None:
                rightids.append(-1) 
            else:
                rightids.append(node.right.id) 
            if node.parent is None:
                parentids.append(-1) 
            else:
                parentids.append(node.parent.id) 
            if node.grid is None:
                gridids.append(-1) 
            else:
                gridids.append(node.grid) 
            if node.split is None:
                splitdims.append(-1)
                splitposs.append(np.nan)
            else:
                splitdims.append(node.split.dim)
                splitposs.append(node.split.pos)

        return nids, parentids, leftids, rightids, les, res, gridids,\
                splitdims, splitposs

    def rebuild_tree_from_array(self, nids, pids, lids,
                               rids, les, res, gids, splitdims, splitposs):
        del self.tree.trunk

        self.tree.trunk = Node(None, 
                    None,
                    None,
                    les[0], res[0], gids[0], nids[0]) 

        N = nids.shape[0]
        for i in xrange(N):
            n = self.get_node(nids[i])
            n.left_edge = les[i]
            n.right_edge = res[i]
            if lids[i] != -1 and n.left is None:
                n.left = Node(n, None, None, None,  
                                      None, None, lids[i])
            if rids[i] != -1 and n.right is None:
                n.right = Node(n, None, None, None, 
                                      None, None, rids[i])
            if gids[i] != -1:
                n.grid = gids[i]

            if splitdims[i] != -1:
                n.split = Split(splitdims[i], splitposs[i])

        mylog.info('AMRKDTree rebuilt, Final Volume: %e' % kd_sum_volume(self.tree.trunk))
        return self.tree.trunk

    def count_volume(self):
        return kd_sum_volume(self.tree.trunk)
    
    def count_cells(self):
        return self.tree.sum_cells() 

if __name__ == "__main__":
    from yt.mods import *
    from time import time
    pf = load('/Users/skillman/simulations/DD1717/DD1717')
    pf.h

    t1 = time()
    hv = AMRKDTree(pf)
    t2 = time()

    print kd_sum_volume(hv.tree.trunk)
    print kd_node_check(hv.tree.trunk)
    print 'Time: %e seconds' % (t2-t1)



