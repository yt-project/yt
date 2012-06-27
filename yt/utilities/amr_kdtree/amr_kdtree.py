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
import numpy as na
from yt.utilities.lib import kdtree_get_choices
from amr_kdtools import kd_is_leaf, kd_sum_volume, kd_node_check, \
        depth_traverse, viewpoint_traverse
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import ParallelAnalysisInterface 
from yt.visualization.volume_rendering.grid_partitioner import HomogenizedVolume
from yt.utilities.lib.grid_traversal import PartitionedGrid
import pdb

def my_break():
    my_debug = False 
    if my_debug: pdb.set_trace()

class Split(object):
    dim = None
    pos = None
    def __init__(self, dim, pos):
        self.dim = dim
        self.pos = pos

class Node(object):
    left = None
    right = None
    parent = None
    data = None
    split = None
    data = None
    def __init__(self, parent, left, right, 
            left_edge, right_edge, grid_id):
        self.left = left
        self.right = right
        self.left_edge = left_edge
        self.right_edge = right_edge
        self.grid = grid_id
        self.parent = parent

class Tree(object):
    trunk = None
    pf = None
    _id_offset = None
    def __init__(self, pf, left=None, right=None):
        self.pf = pf
        self._id_offset = self.pf.h.grids[0]._id_offset
        if left is None:
            left = na.array([-na.inf]*3)
        if right is None:
            right = na.array([na.inf]*3)
        self.trunk = Node(None, None, None,
                left, right, None)
        self.build()

    def build(self, grids = None):
        if grids is None:
            level_iter = self.pf.hierarchy.get_levels()
            while True:
                try:
                    grids = level_iter.next()
                except:
                    break
                gles = na.array([g.LeftEdge for g in grids])
                gres = na.array([g.RightEdge for g in grids])
                gids = na.array([g.id for g in grids])

                add_grids(self.trunk, gles, gres, gids)
                #print 'Checking level %i' % grids[0].Level
                #self.check_tree()

    def check_tree(self):
        for node in depth_traverse(self):
            if node.grid is None: 
                print 'Node has no grid', node, node.left_edge, node.right_edge
                continue
            grid = self.pf.h.grids[node.grid - self._id_offset]
            dds = grid.dds
            gle = grid.LeftEdge
            gre = grid.RightEdge
            li = na.rint((node.left_edge-gle)/dds).astype('int32')
            ri = na.rint((node.right_edge-gle)/dds).astype('int32')
            dims = (ri - li).astype('int32')
            assert(na.all(grid.LeftEdge <= node.left_edge))
            assert(na.all(grid.RightEdge >= node.right_edge))
            print grid, dims, li, ri

    def insert_level(self, level, node):
        if node.left is not None:
            self.insert_level(level, node.left)
        if node.right is not None:
            self.insert_level(level, node.right)
        if node.grid is None: return
        if self.pf.h.grid_levels[node.grid-self._id_offset][0] == (level-1):
            self.insert_children(node)

    def insert_children(self, node):
        thisgrid = self.pf.h.grids[node.grid - self._id_offset]
        if len(thisgrid.Children) == 0: return
        children = [child for child in thisgrid.Children
                    if na.all(child.LeftEdge < node.right_edge) &
                    na.all(child.RightEdge > node.left_edge)]

        # If we have children, get all the new grids, and keep building the tree
        if len(children) > 0:
            gles = na.array([g.LeftEdge for g in children])
            gres = na.array([g.RightEdge for g in children])
            gids = na.array([g.id for g in children])
            insert_grids(node, gles, gres, gids)

def add_grids(node, gles, gres, gids):
    if kd_is_leaf(node):
        insert_grids(node, gles, gres, gids)
    else:
        less_ids = gles[:,node.split.dim] < node.split.pos
        if len(less_ids) > 0:
            add_grids(node.left, gles[less_ids], gres[less_ids], gids[less_ids])

        greater_ids = gres[:,node.split.dim] > node.split.pos
        if len(greater_ids) > 0:
            add_grids(node.right, gles[greater_ids], gres[greater_ids], gids[greater_ids])

def insert_grids(node, gles, gres, grid_ids):
    if len(grid_ids) < 1:
        #node.grid = node.parent.grid
        assert(node.grid is not None)
        return

    if len(grid_ids) == 1:
        if na.all(gles[0] <= node.left_edge) and \
                na.all(gres[0] >= node.right_edge):
            node.grid = grid_ids[0]
            #print 'Created a node using grid %i with le, re' % node.grid, node.left_edge, node.right_edge
            assert(node.grid is not None)
            return

    # Find a Split
    data = na.array([(gles[i,:], gres[i,:]) for i in xrange(grid_ids.shape[0])], copy=False)
    best_dim, split_pos, less_ids, greater_ids = \
        kdtree_get_choices(data, node.left_edge, node.right_edge)
    split = Split(best_dim, split_pos)

    del data, best_dim, split_pos

    # Create a Split
    divide(node, split)

    # Populate Left Node
    #print 'Inserting left node', node.left_edge, node.right_edge
    insert_grids(node.left, gles[less_ids], gres[less_ids], grid_ids[less_ids])

    # Populate Right Node
    #print 'Inserting right node', node.left_edge, node.right_edge
    insert_grids(node.right, gles[greater_ids], gres[greater_ids], grid_ids[greater_ids])

    del less_ids, greater_ids
    return

class AMRKDTree(HomogenizedVolume):
    current_vcds = []
    tree = None
    no_ghost = True
    fields = ['Density']
    log_fields = [True]
    _id_offset = None
    current_saved_grids = []
    pf = None
    bricks = []
    brick_dimensions = [] 
    _initialized = False
    def __init__(self, pf,  l_max=None, le=None, re=None,
                 fields=None, no_ghost=False,
                 tree_type='domain',log_fields=None, merge_trees=False):

        ParallelAnalysisInterface.__init__(self)

        self.pf = pf
        self.l_max = l_max
        if fields is None: fields = ["Density"]
        self.fields = ensure_list(fields)

        self.no_ghost = no_ghost
        if log_fields is not None:
            log_fields = ensure_list(log_fields)
        else:
            log_fields = [self.pf.field_info[field].take_log
                         for field in self.fields]

        self.log_fields = log_fields
        self._id_offset = pf.h.grids[0]._id_offset

        if le is None:
            self.le = pf.domain_left_edge
        else:
            self.le = na.array(le)
        if re is None:
            self.re = pf.domain_right_edge
        else:
            self.re = na.array(re)

        print 'Building tree with le,re ', self.le, self.re
        self.tree = Tree(pf, self.le, self.re)

    def initialize_source(self):
        if self._initialized : return
        bricks = []
        for b in self.traverse():
            bricks.append(b)
        self.bricks = na.array(bricks)
        self.brick_dimensions = na.array(self.brick_dimensions)
        self._initialized = True

    def traverse(self, viewpoint=None):
        if viewpoint is None:
            for node in depth_traverse(self.tree):
                if kd_is_leaf(node):
                    yield self.get_brick_data(node)
        else:
            for node in viewpoint_traverse(self.tree, viewpoint):
                if kd_is_leaf(node):
                    yield self.get_brick_data(node)

    def get_brick_data(self, node):
        if node.data is not None: return node.data
        grid = self.pf.h.grids[node.grid - self._id_offset]
        dds = grid.dds
        gle = grid.LeftEdge
        gre = grid.RightEdge
        li = na.rint((node.left_edge-gle)/dds).astype('int32')
        ri = na.rint((node.right_edge-gle)/dds).astype('int32')
        dims = (ri - li).astype('int32')
        assert(na.all(grid.LeftEdge <= node.left_edge))
        assert(na.all(grid.RightEdge >= node.right_edge))

        if grid in self.current_saved_grids:
            dds = self.current_vcds[self.current_saved_grids.index(grid)]
        else:
            dds = []
            for i,field in enumerate(self.fields):
                vcd = grid.get_vertex_centered_data(field,smoothed=True,no_ghost=self.no_ghost).astype('float64')
                if self.log_fields[i]: vcd = na.log10(vcd)
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
        if not self._initialized: self.brick_dimensions.append(dims)
        return brick

def new_right(Node, split):
    new_right = na.empty(3, dtype='float64')
    new_right[:] = Node.right_edge[:]
    new_right[split.dim] = split.pos
    return new_right

def new_left(Node, split):
    new_left = na.empty(3, dtype='float64')
    new_left[:] = Node.left_edge[:]
    new_left[split.dim] = split.pos
    return new_left

def divide(node, split):
    # Create a Split
    node.split = split
    node.left = Node(node, None, None,
            node.left_edge, new_right(node, split), node.grid)
    node.right = Node(node, None, None,
            new_left(node, split), node.right_edge, node.grid)
    return


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


