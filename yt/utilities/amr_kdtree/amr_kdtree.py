"""
AMR kD-Tree Framework


"""
from __future__ import print_function
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import operator
import numpy as np

from yt.funcs import \
    iterable, \
    mylog
from yt.utilities.on_demand_imports import _h5py as h5py
from yt.utilities.amr_kdtree.amr_kdtools import \
    receive_and_reduce, \
    send_to_parent, \
    scatter_image
from yt.utilities.lib.amr_kdtools import Node
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.lib.partitioned_grid import PartitionedGrid
from yt.utilities.math_utils import periodic_position
from yt.geometry.grid_geometry_handler import GridIndex

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

def _apply_log(data, log_changed, log_new):
    '''Helper used to set log10/10^ to data in AMRKDTree'''
    if not log_changed:
        return
    if log_new:
        np.log10(data, data)
    else:
        np.power(10.0, data, data)

class Tree(object):
    def __init__(self, ds, comm_rank=0, comm_size=1, left=None, right=None,
        min_level=None, max_level=None, data_source=None):

        self.ds = ds
        try:
            self._id_offset = ds.index.grids[0]._id_offset
        except AttributeError:
            self._id_offset = 0

        if data_source is None:
            data_source = ds.all_data()
        self.data_source = data_source
        if left is None:
            left = np.array([-np.inf]*3)
        if right is None:
            right = np.array([np.inf]*3)

        if min_level is None: min_level = 0
        if max_level is None: max_level = ds.index.max_level
        self.min_level = min_level
        self.max_level = max_level
        self.comm_rank = comm_rank
        self.comm_size = comm_size
        self.trunk = Node(None, None, None, left, right, -1, 1)
        self.build()

    def add_grids(self, grids):
        gles = np.array([g.LeftEdge for g in grids])
        gres = np.array([g.RightEdge for g in grids])
        gids = np.array([g.id for g in grids], dtype="int64")
        self.trunk.add_grids(gids.size, gles, gres, gids,
                    self.comm_rank, self.comm_size)
        del gles, gres, gids, grids

    def build(self):
        lvl_range = range(self.min_level, self.max_level+1)
        for lvl in lvl_range:
            #grids = self.data_source.select_grids(lvl)
            grids = np.array([b for b, mask in self.data_source.blocks if b.Level == lvl])
            if len(grids) == 0: continue
            self.add_grids(grids)

    def check_tree(self):
        for node in self.trunk.depth_traverse():
            if node.grid == -1:
                continue
            grid = self.ds.index.grids[node.grid - self._id_offset]
            dds = grid.dds
            gle = grid.LeftEdge
            nle = self.ds.arr(node.get_left_edge(), input_units="code_length")
            nre = self.ds.arr(node.get_right_edge(), input_units="code_length")
            li = np.rint((nle-gle)/dds).astype('int32')
            ri = np.rint((nre-gle)/dds).astype('int32')
            dims = (ri - li).astype('int32')
            assert(np.all(grid.LeftEdge <= nle))
            assert(np.all(grid.RightEdge >= nre))
            assert(np.all(dims > 0))
            # print grid, dims, li, ri

        # Calculate the Volume
        vol = self.trunk.kd_sum_volume()
        mylog.debug('AMRKDTree volume = %e' % vol)
        self.trunk.kd_node_check()

    def sum_cells(self, all_cells=False):
        cells = 0
        for node in self.trunk.depth_traverse():
            if node.grid == -1:
                continue
            if not all_cells and not node.kd_is_leaf():
                continue
            grid = self.ds.index.grids[node.grid - self._id_offset]
            dds = grid.dds
            gle = grid.LeftEdge
            nle = self.ds.arr(node.get_left_edge(), input_units="code_length")
            nre = self.ds.arr(node.get_right_edge(), input_units="code_length")
            li = np.rint((nle-gle)/dds).astype('int32')
            ri = np.rint((nre-gle)/dds).astype('int32')
            dims = (ri - li).astype('int32')
            cells += np.prod(dims)
        return cells


class AMRKDTree(ParallelAnalysisInterface):
    r"""A KDTree for AMR data. 

    Not applicable to particle or octree-based datasets.

    """

    fields = None
    log_fields = None
    no_ghost = True

    def __init__(self, ds, min_level=None, max_level=None,
                 data_source=None):

        if not issubclass(ds.index.__class__, GridIndex):
            raise RuntimeError("AMRKDTree does not support particle or octree-based data.")

        ParallelAnalysisInterface.__init__(self)

        self.ds = ds
        self.current_vcds = []
        self.current_saved_grids = []
        self.bricks = []
        self.brick_dimensions = []
        self.sdx = ds.index.get_smallest_dx()

        self._initialized = False
        try:
            self._id_offset = ds.index.grids[0]._id_offset
        except AttributeError:
            self._id_offset = 0

        if data_source is None:
            data_source = self.ds.all_data()
        self.data_source = data_source

        mylog.debug('Building AMRKDTree')
        self.tree = Tree(ds, self.comm.rank, self.comm.size,
                         min_level=min_level, max_level=max_level,
                         data_source=data_source)

    def set_fields(self, fields, log_fields, no_ghost, force=False):
        new_fields = self.data_source._determine_fields(fields)
        regenerate_data = self.fields is None or \
                          len(self.fields) != len(new_fields) or \
                          self.fields != new_fields or force
        if not iterable(log_fields):
            log_fields = [log_fields]
        new_log_fields = list(log_fields)
        self.tree.trunk.set_dirty(regenerate_data)
        self.fields = new_fields

        if self.log_fields is not None and not regenerate_data:
            flip_log = list(map(operator.ne, self.log_fields, new_log_fields))
        else:
            flip_log = [False] * len(new_log_fields)
        self.log_fields = new_log_fields

        self.no_ghost = no_ghost
        del self.bricks, self.brick_dimensions
        self.brick_dimensions = []
        bricks = []

        for b in self.traverse():
            list(map(_apply_log, b.my_data, flip_log, self.log_fields))
            bricks.append(b)
        self.bricks = np.array(bricks)
        self.brick_dimensions = np.array(self.brick_dimensions)
        self._initialized = True

    def initialize_source(self, fields, log_fields, no_ghost):
        if fields == self.fields and log_fields == self.log_fields and \
                no_ghost == self.no_ghost:
            return
        self.set_fields(fields, log_fields, no_ghost)

    def traverse(self, viewpoint=None):
        for node in self.tree.trunk.kd_traverse(viewpoint=viewpoint):
            yield self.get_brick_data(node)

    def slice_traverse(self, viewpoint = None):
        if not hasattr(self.ds.index, "grid"):
            raise NotImplementedError
        for node in self.tree.trunk.kd_traverse(viewpoint=viewpoint):
            grid = self.ds.index.grids[node.grid - self._id_offset]
            dds = grid.dds
            gle = grid.LeftEdge.in_units("code_length").ndarray_view()
            nle = node.get_left_edge()
            nre = node.get_right_edge()
            li = np.rint((nle-gle)/dds).astype('int32')
            ri = np.rint((nre-gle)/dds).astype('int32')
            dims = (ri - li).astype('int32')
            sl = (slice(li[0], ri[0]),
                  slice(li[1], ri[1]),
                  slice(li[2], ri[2]))
            gi = grid.get_global_startindex() + li
            yield grid, node, (sl, dims, gi)

    def get_node(self, nodeid):
        path = np.binary_repr(nodeid)
        depth = 1
        temp = self.tree.trunk
        for depth in range(1, len(path)):
            if path[depth] == '0':
                temp = temp.left
            else:
                temp = temp.right
        assert(temp is not None)
        return temp

    def locate_node(self, pos):
        return self.tree.trunk.find_node(pos)

    def get_reduce_owners(self):
        owners = {}
        for bottom_id in range(self.comm.size, 2*self.comm.size):
            temp = self.get_node(bottom_id)
            owners[temp.node_id] = temp.node_id - self.comm.size
            while temp is not None:
                if temp.parent is None: break
                if temp == temp.parent.right:
                    break
                temp = temp.parent
                owners[temp.node_id] = owners[temp.left.node_id]
        return owners

    def reduce_tree_images(self, image, viewpoint):
        if self.comm.size <= 1: return image
        myrank = self.comm.rank
        nprocs = self.comm.size
        owners = self.get_reduce_owners()
        node = self.get_node(nprocs + myrank)

        while owners[node.parent.node_id] == myrank:
            split_dim = node.parent.get_split_dim()
            split_pos = node.parent.get_split_pos()
            add_to_front = viewpoint[split_dim] >= split_pos
            image = receive_and_reduce(self.comm,
                                       owners[node.parent.right.node_id],
                                       image, add_to_front)
            if node.parent.node_id == 1: break
            else: node = node.parent
        else:
            send_to_parent(self.comm, owners[node.parent.node_id], image)

        return scatter_image(self.comm, owners[1], image)

    def get_brick_data(self, node):
        if node.data is not None and not node.dirty:
            return node.data
        grid = self.ds.index.grids[node.grid - self._id_offset]
        dds = grid.dds.ndarray_view()
        gle = grid.LeftEdge.ndarray_view()
        nle = node.get_left_edge()
        nre = node.get_right_edge()
        li = np.rint((nle-gle)/dds).astype('int32')
        ri = np.rint((nre-gle)/dds).astype('int32')
        dims = (ri - li).astype('int32')
        assert(np.all(grid.LeftEdge <= nle))
        assert(np.all(grid.RightEdge >= nre))

        if grid in self.current_saved_grids and not node.dirty:
            dds = self.current_vcds[self.current_saved_grids.index(grid)]
        else:
            dds = []
            vcd = grid.get_vertex_centered_data(self.fields, smoothed=True,
                                                no_ghost=self.no_ghost)
            for i, field in enumerate(self.fields):
                if self.log_fields[i]:
                    dds.append(np.log10(vcd[field].astype('float64')))
                else:
                    dds.append(vcd[field].astype('float64'))
                self.current_saved_grids.append(grid)
                self.current_vcds.append(dds)

        if self.data_source.selector is None:
            mask = np.ones(dims, dtype='uint8')
        else:
            mask = self.data_source.selector.fill_mask(grid)[li[0]:ri[0], li[1]:ri[1], li[2]:ri[2] ].astype('uint8')

        data = [d[li[0]:ri[0]+1,
                  li[1]:ri[1]+1,
                  li[2]:ri[2]+1].copy() for d in dds]

        brick = PartitionedGrid(grid.id, data,
                                mask,
                                nle.copy(),
                                nre.copy(),
                                dims.astype('int64'))
        node.data = brick
        node.dirty = False
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
        new_positions = [periodic_position(p, self.ds) for p in new_positions]
        grids[in_grid] = grid

        get_them = np.argwhere(in_grid).ravel()
        cis[in_grid] = new_cis[in_grid]

        if (in_grid).sum()>0:
            grids[np.logical_not(in_grid)] = \
                [self.ds.index.grids[self.locate_brick(new_positions[i]).grid -
                                 self._id_offset]
                 for i in get_them]
            cis[np.logical_not(in_grid)] = \
                [(new_positions[i]-grids[i].LeftEdge)/
                 grids[i].dds for i in get_them]
        cis = [tuple(_ci) for _ci in cis]
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
        grid = self.ds.index.grids[self.locate_brick(position).grid -
                               self._id_offset]
        ci = ((position-grid.LeftEdge)/grid.dds).astype('int64')
        return self.locate_neighbors(grid,ci)

    def store_kd_bricks(self, fn=None):
        if not self._initialized:
            self.initialize_source()
        if fn is None:
            fn = '%s_kd_bricks.h5'%self.ds
        if self.comm.rank != 0:
            self.comm.recv_array(self.comm.rank-1, tag=self.comm.rank-1)
        f = h5py.File(fn,'w')
        for node in self.tree.depth_traverse():
            i = node.node_id
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
            fn = '%s_kd_bricks.h5' % self.ds
        if self.comm.rank != 0:
            self.comm.recv_array(self.comm.rank-1, tag=self.comm.rank-1)
        try:
            f = h5py.File(fn,"a")
            for node in self.tree.depth_traverse():
                i = node.node_id
                if node.grid != -1:
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
        self.rebuild_tree_from_array(nid, pid, lid,
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
        for node in self.tree.trunk.depth_first_touch():
            nids.append(node.node_id)
            les.append(node.get_left_edge())
            res.append(node.get_right_edge())
            if node.left is None:
                leftids.append(-1)
            else:
                leftids.append(node.left.node_id)
            if node.right is None:
                rightids.append(-1)
            else:
                rightids.append(node.right.node_id)
            if node.parent is None:
                parentids.append(-1)
            else:
                parentids.append(node.parent.node_id)
            if node.grid is None:
                gridids.append(-1)
            else:
                gridids.append(node.grid)
            splitdims.append(node.get_split_dim())
            splitposs.append(node.get_split_pos())

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
        for i in range(N):
            n = self.get_node(nids[i])
            n.set_left_edge(les[i])
            n.set_right_edge(res[i])
            if lids[i] != -1 and n.left is None:
                n.left = Node(n, None, None,
                              np.zeros(3, dtype='float64'),
                              np.zeros(3, dtype='float64'),
                              -1, lids[i])
            if rids[i] != -1 and n.right is None:
                n.right = Node(n, None, None,
                               np.zeros(3, dtype='float64'),
                               np.zeros(3, dtype='float64'),
                               -1, rids[i])
            if gids[i] != -1:
                n.grid = gids[i]

            if splitdims[i] != -1:
                n.create_split(splitdims[i], splitposs[i])

        mylog.info('AMRKDTree rebuilt, Final Volume: %e' % self.tree.trunk.kd_sum_volume())
        return self.tree.trunk

    def count_volume(self):
        return self.tree.trunk.kd_sum_volume()

    def count_cells(self):
        return self.tree.sum_cells()

if __name__ == "__main__":
    import yt
    from time import time
    ds = yt.load('/Users/skillman/simulations/DD1717/DD1717')
    ds.index

    t1 = time()
    hv = AMRKDTree(ds)
    t2 = time()

    print(hv.tree.trunk.kd_sum_volume())
    print(hv.tree.trunk.kd_node_check())
    print('Time: %e seconds' % (t2-t1))
