"""
AMR kD-Tree Tools 

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
import numpy as np
from yt.funcs import *
from yt.utilities.lib import kdtree_get_choices

def _lchild_id(node_id): return (node_id<<1)
def _rchild_id(node_id): return (node_id<<1) + 1
def _parent_id(node_id): return (node_id-1) >> 1

class Node(object):
    def __init__(self, parent, left, right,
            left_edge, right_edge, grid_id, node_id):
        self.left = left
        self.right = right
        self.left_edge = left_edge
        self.right_edge = right_edge
        self.grid = grid_id
        self.parent = parent
        self.id = node_id
        self.data = None

class Split(object):
    def __init__(self, dim, pos):
        self.dim = dim
        self.pos = pos

def should_i_build(node, rank, size):
    if (node.id < size) or (node.id >= 2*size):
        return True
    elif node.id - size == rank:
        return True
    else:
        return False

def add_grids(node, gles, gres, gids, rank, size):
    if not should_i_build(node, rank, size):
        return

    if kd_is_leaf(node):
        insert_grids(node, gles, gres, gids, rank, size)
    else:
        less_ids = gles[:,node.split.dim] < node.split.pos
        if len(less_ids) > 0:
            add_grids(node.left, gles[less_ids], gres[less_ids],
                      gids[less_ids], rank, size)

        greater_ids = gres[:,node.split.dim] > node.split.pos
        if len(greater_ids) > 0:
            add_grids(node.right, gles[greater_ids], gres[greater_ids],
                      gids[greater_ids], rank, size)

def should_i_split(node, rank, size):
    return node.id < size

def geo_split(node, gles, gres, grid_ids, rank, size):
    big_dim = np.argmax(gres[0]-gles[0])
    new_pos = (gres[0][big_dim] + gles[0][big_dim])/2.
    old_gre = gres[0].copy()
    new_gle = gles[0].copy()
    new_gle[big_dim] = new_pos
    gres[0][big_dim] = new_pos
    gles = np.append(gles, np.array([new_gle]), axis=0)
    gres = np.append(gres, np.array([old_gre]), axis=0)
    grid_ids = np.append(grid_ids, grid_ids, axis=0)

    split = Split(big_dim, new_pos)

    # Create a Split
    divide(node, split)

    # Populate Left Node
    #print 'Inserting left node', node.left_edge, node.right_edge
    insert_grids(node.left, gles[:1], gres[:1],
            grid_ids[:1], rank, size)

    # Populate Right Node
    #print 'Inserting right node', node.left_edge, node.right_edge
    insert_grids(node.right, gles[1:], gres[1:],
            grid_ids[1:], rank, size)
    return

def insert_grids(node, gles, gres, grid_ids, rank, size):
    if not should_i_build(node, rank, size) or grid_ids.size == 0:
        return

    if len(grid_ids) == 1:
        # If we should continue to split based on parallelism, do so!
        if should_i_split(node, rank, size):
            geo_split(node, gles, gres, grid_ids, rank, size)
            return

        if np.all(gles[0] <= node.left_edge) and \
                np.all(gres[0] >= node.right_edge):
            node.grid = grid_ids[0]
            assert(node.grid is not None)
            return

    # Split the grids
    check = split_grids(node, gles, gres, grid_ids, rank, size)
    # If check is -1, then we have found a place where there are no choices.
    # Exit out and set the node to None.
    if check == -1:
        node.grid = None
    return

def split_grids(node, gles, gres, grid_ids, rank, size):
    # Find a Split
    data = np.array([(gles[i,:], gres[i,:]) for i in
        xrange(grid_ids.shape[0])], copy=False)
    best_dim, split_pos, less_ids, greater_ids = \
        kdtree_get_choices(data, node.left_edge, node.right_edge)

    # If best_dim is -1, then we have found a place where there are no choices.
    # Exit out and set the node to None.
    if best_dim == -1:
        return -1

    split = Split(best_dim, split_pos)

    del data, best_dim, split_pos

    # Create a Split
    divide(node, split)

    # Populate Left Node
    #print 'Inserting left node', node.left_edge, node.right_edge
    insert_grids(node.left, gles[less_ids], gres[less_ids],
                 grid_ids[less_ids], rank, size)

    # Populate Right Node
    #print 'Inserting right node', node.left_edge, node.right_edge
    insert_grids(node.right, gles[greater_ids], gres[greater_ids],
                 grid_ids[greater_ids], rank, size)

    return

def new_right(Node, split):
    new_right = Node.right_edge.copy()
    new_right[split.dim] = split.pos
    return new_right

def new_left(Node, split):
    new_left = Node.left_edge.copy()
    new_left[split.dim] = split.pos
    return new_left

def divide(node, split):
    # Create a Split
    node.split = split
    node.left = Node(node, None, None,
            node.left_edge, new_right(node, split), node.grid,
                     _lchild_id(node.id))
    node.right = Node(node, None, None,
            new_left(node, split), node.right_edge, node.grid,
                      _rchild_id(node.id))
    return

def kd_sum_volume(node):
    if (node.left is None) and (node.right is None):
        if node.grid is None:
            return 0.0
        return np.prod(node.right_edge - node.left_edge)
    else:
        return kd_sum_volume(node.left) + kd_sum_volume(node.right)

def kd_node_check(node):
    assert (node.left is None) == (node.right is None)
    if (node.left is None) and (node.right is None):
        if node.grid is not None:
            return np.prod(node.right_edge - node.left_edge)
        else: return 0.0
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

def depth_traverse(tree, max_node=None):
    '''
    Yields a depth-first traversal of the kd tree always going to
    the left child before the right.
    '''
    current = tree.trunk
    previous = None
    if max_node is None:
        max_node = np.inf
    while current is not None:
        yield current
        current, previous = step_depth(current, previous)
        if current is None: break
        if current.id >= max_node:
            current = current.parent
            previous = current.right

def depth_first_touch(tree, max_node=None):
    '''
    Yields a depth-first traversal of the kd tree always going to
    the left child before the right.
    '''
    current = tree.trunk
    previous = None
    if max_node is None:
        max_node = np.inf
    while current is not None:
        if previous is None or previous.parent != current:
            yield current
        current, previous = step_depth(current, previous)
        if current is None: break
        if current.id >= max_node:
            current = current.parent
            previous = current.right

def breadth_traverse(tree):
    '''
    Yields a breadth-first traversal of the kd tree always going to
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


def receive_and_reduce(comm, incoming_rank, image, add_to_front):
    mylog.debug( 'Receiving image from %04i' % incoming_rank)
    #mylog.debug( '%04i receiving image from %04i'%(self.comm.rank,back.owner))
    arr2 = comm.recv_array(incoming_rank, incoming_rank).reshape(
        (image.shape[0], image.shape[1], image.shape[2]))

    if add_to_front:
        front = arr2
        back = image
    else:
        front = image
        back = arr2

    if image.shape[2] == 3:
        # Assume Projection Camera, Add
        np.add(image, front, image)
        return image

    ta = 1.0 - front[:,:,3]
    np.maximum(ta, 0.0, ta)
    # This now does the following calculation, but in a memory
    # conservative fashion
    # image[:,:,i  ] = front[:,:,i] + ta*back[:,:,i]
    image = back.copy()
    for i in range(4):
        np.multiply(image[:,:,i], ta, image[:,:,i])
    np.add(image, front, image)
    return image

def send_to_parent(comm, outgoing_rank, image):
    mylog.debug( 'Sending image to %04i' % outgoing_rank)
    comm.send_array(image, outgoing_rank, tag=comm.rank)

def scatter_image(comm, root, image):
    mylog.debug( 'Scattering from %04i' % root)
    image = comm.mpi_bcast(image, root=root)
    return image

def find_node(node, pos):
    """
    Find the AMRKDTree node enclosing a position
    """
    assert(np.all(node.left_edge <= pos))
    assert(np.all(node.right_edge > pos))
    while not kd_is_leaf(node):
        if pos[node.split.dim] < node.split.pos:
            node = node.left
        else:
            node = node.right
    return node

