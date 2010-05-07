"""
AMR kD-Tree Framework

Authors: Samuel Skillman <samskillman@gmail.com>
         Wil St. Charles <fallen751@gmail.com>
Affiliation: University of Colorado at Boulder
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Samuel Skillman.  All Rights Reserved.

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
from copy import * 
from yt.mods import *
from time import time

class AMRKDTree(object):
    def __init__(self, pf,  L_MAX=None, extra_dim_cost=0, le=None, re=None):
        self.pf = pf
        if L_MAX is None:
            self.L_MAX = self.pf.hierarchy.max_level+1
        else:
            self.L_MAX = np.min([L_MAX,self.pf.hierarchy.max_level+1])

        self.extra_dim_cost=extra_dim_cost
        if le is None: self.domain_left_edge = np.array([0.0]*3)
        else: self.domain_left_edge = le
        if re is None: self.domain_right_edge = np.array([1.0]*3)
        else: self.domain_right_edge = re
        print 'Making kd tree from le ', self.domain_left_edge, 'to ', self.domain_right_edge

        self.root_grids = pf.hierarchy.get_levels().next()
        self.leaf_count = 0
        self.current_parent = -1
        self.dim = 0
        tt1 = time()
        self.tree = self.__build(self.root_grids, None, self.domain_left_edge, self.domain_right_edge)
        print 'It took %e seconds to build the tree' % (time()-tt1)
        self.total_cost = self.count_cost(self.tree)
        
    class node(object):
        pass

    class leafnode(node):
        def __init__(self,leaf_id, grid_id, leaf_l_corner, leaf_r_corner):
            self.leaf_id = leaf_id
            self.grid = grid_id
            self.leaf_l_corner = leaf_l_corner
            self.leaf_r_corner = leaf_r_corner
            dds = self.grid.dds
            gle = self.grid.LeftEdge
            gre = self.grid.RightEdge
            self.li = ((leaf_l_corner-gle)/dds).astype('int32')
            self.ri = ((leaf_r_corner-gle)/dds).astype('int32')
            self.dims = (self.ri - self.li).astype('int32')
            # Here the cost is actually inversely proportional to 4**Level (empirical)
            self.cost = (np.prod(self.dims)/4.**self.grid.Level).astype('int64')
            # Here is the old way
            # self.cost = np.prod(self.dims).astype('int64')
            self.owner = -1
            del dds, gle, gre
            
    class dividing_node(node):
        def __init__(self, split_ax, split_pos, left_children, right_children):
            self.split_ax = split_ax
            self.split_pos = split_pos
            self.left_children = left_children
            self.right_children = right_children
            self.cost = 0.0
            self.owner = -1

    def count_cost(self,node):
        if isinstance(node,AMRKDTree.leafnode):
            return node.cost
        else:
            node.cost = self.count_cost(node.left_children) + self.count_cost(node.right_children)
            return node.cost

    def __build(self, grids, parent, l_corner, r_corner):
        if len(grids) == 0:
            self.leaf_count += 1
            return AMRKDTree.leafnode(self.leaf_count-1, parent, l_corner, r_corner)

        if len(grids) == 1:
            thisgrid = grids[0]
            if  ( (thisgrid.LeftEdge[0] <= l_corner[0]) and (thisgrid.RightEdge[0] >= r_corner[0]) and
                  (thisgrid.LeftEdge[1] <= l_corner[1]) and (thisgrid.RightEdge[1] >= r_corner[1]) and
                  (thisgrid.LeftEdge[2] <= l_corner[2]) and (thisgrid.RightEdge[2] >= r_corner[2]) ):
                children = []
                if (len(thisgrid.Children) > 0) and (thisgrid.Level < self.L_MAX):
                    children = np.unique(thisgrid._get_child_index_mask()[thisgrid._get_child_index_mask() >= 0])
                    children = self.pf.hierarchy.grids[children -1]
                    child_l_data = np.array([child.LeftEdge for child in children])
                    child_r_data = np.array([child.RightEdge for child in children])
                    children_we_want =  (child_l_data[:,0] < r_corner[0])*(child_r_data [:,0] > l_corner[0])* \
                        (child_l_data[:,1] < r_corner[1])*(child_r_data [:,1] > l_corner[1])* \
                        (child_l_data[:,2] < r_corner[2])*(child_r_data [:,2] > l_corner[2])
                    children = children[children_we_want]
                return self.__build(children, thisgrid, l_corner, r_corner)

        l_data = np.array([child.LeftEdge for child in grids])
        r_data = np.array([child.RightEdge for child in grids])

        best_dim = 0
        best_choices = []
        for d in range(3):
            choices = np.unique(np.clip(np.concatenate(([l_corner[d]],l_data[:,d],r_data[:,d],[r_corner[d]])),l_corner[d],r_corner[d]))
            if len(choices) > len(best_choices):
                best_choices = choices
                best_dim = d
        best_choices.sort()
        split = best_choices[len(best_choices)/2]

        less_ids = np.nonzero(l_data[:,best_dim] < split)
        greater_ids = np.nonzero(split < r_data[:,best_dim])
        return AMRKDTree.dividing_node(best_dim, split, 
                                       self.__build(grids[less_ids], parent,
                                                    l_corner, self.__corner_bounds(best_dim, split, current_right=r_corner)),
                                       self.__build(grids[greater_ids], parent,
                                                    self.__corner_bounds(best_dim, split, current_left=l_corner), r_corner))
        
    def __corner_bounds(self, split_dim, split, current_left = None, current_right = None):
        if(current_left is not None):
            new_left = deepcopy(current_left)
            new_left[split_dim] = split
            return new_left
        elif(current_right is not None):
            new_right = deepcopy(current_right)
            new_right[split_dim] = split
            return new_right
