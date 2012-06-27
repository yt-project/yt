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
import numpy as na
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
        node.grid = node.parent.grid
        assert(node.grid is not None)
        return

    if len(grid_ids) == 1:
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


def kd_sum_volume(node):
    if (node.left is None)  and (node.right is None):
        return na.prod(node.right_edge - node.left_edge)
    else:
        return kd_sum_volume(node.left) + kd_sum_volume(node.right)

def kd_node_check(node):
    assert (node.left is None) == (node.right is None)
    if (node.left is None)  and (node.right is None):
        assert(node.grid is not None)
        return na.prod(node.right_edge - node.left_edge)
    else:
        return kd_node_check(node.left)+kd_node_check(node.right)

def kd_is_leaf(node):
    has_l_child = node.left is None
    has_r_child = node.right is None
    assert has_l_child == has_r_child
    return has_l_child

def step_depth(current, previous):
    '''
    Takes a single step in the depth-first traversal
    '''
    if kd_is_leaf(current): # At a leaf, move back up
        previous = current
        current = current.parent

    elif current.parent is previous: # Moving down, go left first
        previous = current
        if current.left is not None:
            current = current.left
        elif current.right is not None:
            current = current.right
        else:
            current = current.parent

    elif current.left is previous: # Moving up from left, go right 
        previous = current
        if current.right is not None:
            current = current.right
        else:
            current = current.parent

    elif current.right is previous: # Moving up from right child, move up
        previous = current
        current = current.parent

    return current, previous

def depth_traverse(tree):
    '''
    Yields a depth-first traversal of the kd tree always going to
    the left child before the right.
    '''
    current = tree.trunk
    previous = None
    while current is not None:
        yield current
        current, previous = step_depth(current, previous)

def viewpoint_traverse(tree, viewpoint):
    '''
    Yields a viewpoint dependent traversal of the kd-tree.  Starts
    with nodes furthest away from viewpoint.
    '''

    current = tree.trunk
    previous = None
    while current is not None:
        yield current
        current, previous = step_viewpoint(current, previous, viewpoint)

def step_viewpoint(current, previous, viewpoint):
    '''
    Takes a single step in the viewpoint based traversal.  Always
    goes to the node furthest away from viewpoint first.
    '''
    if kd_is_leaf(current): # At a leaf, move back up
        previous = current
        current = current.parent
    elif current.split.dim is None: # This is a dead node
        previous = current
        current = current.parent

    elif current.parent is previous: # Moving down
        previous = current
        if viewpoint[current.split.dim] <= current.split.pos:
            if current.right is not None:
                current = current.right
            else:
                previous = current.right
        else:
            if current.left is not None:
                current = current.left
            else:
                previous = current.left

    elif current.right is previous: # Moving up from right 
        previous = current
        if viewpoint[current.split.dim] <= current.split.pos:
            if current.left is not None:
                current = current.left
            else:
                current = current.parent
        else:
            current = current.parent

    elif current.left is previous: # Moving up from left child
        previous = current
        if viewpoint[current.split.dim] > current.split.pos:
            if current.right is not None:
                current = current.right
            else:
                current = current.parent
        else:
            current = current.parent

    return current, previous

