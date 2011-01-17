"""
AMR kD-Tree Framework

Authors: Samuel Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
Wil St. Charles <fallen751@gmail.com>
Affiliation: University of Colorado at Boulder

Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Samuel Skillman.  All Rights Reserved.

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
from yt.visualization.volume_rendering.grid_partitioner import HomogenizedVolume
from yt.utilities.amr_utils import PartitionedGrid
from yt.utilities.performance_counters import yt_counters, time_function
import yt.utilities.parallel_tools.parallel_analysis_interface as PT

from yt.config import ytcfg

from time import time
import h5py
my_rank = ytcfg.getint("yt", "__parallel_rank")
nprocs = ytcfg.getint("yt", "__parallel_size")

def corner_bounds(split_dim, split, current_left = None, current_right = None):
    r"""
    Given a kd-Tree split dimension and position and bound to be
    modified, returns the new bound.

    A simple function that replaces the `split_dim` dimension of the
    current left or right bound with `split`. Left or Right bound is
    chosen by specifying the `current_left` or `current_right`.
    """
    if(current_left is not None):
        new_left = na.array([current_left[0],current_left[1],current_left[2]])
        new_left[split_dim] = split
        return new_left
    elif(current_right is not None):
        new_right = na.array([current_right[0],current_right[1],current_right[2]])
        new_right[split_dim] = split
        return new_right

def _lchild_id(id): return (id<<1) + 1
def _rchild_id(id): return (id<<1) + 2
def _parent_id(id): return (id-1)>>1

class MasterNode(object):
    r"""
    A MasterNode object is the building block of the AMR kd-Tree.
    Used during the construction to act as both dividing nodes and
    leaf nodes.
    """
    def __init__(self):
        self.grids = None
        self.parent = None
        self.parent_grid = None
        self.l_corner = None
        self.r_corner = None

        self.split_ax = None
        self.split_pos = None

        self.left_children = None
        self.right_children = None
        self.cost = 0
        self.owner = -1

        self.id = None
        
        self.grid = None
        self.brick = None
        self.li = None
        self.ri = None
        self.dims = None

        self.done = 0
        self.cast_done = 0

def set_leaf(thisnode, grid_id, leaf_l_corner, leaf_r_corner):
    r"""
    Sets leaf properties.

    Parameters
    ----------
    thisnode : `MasterNode`
        AMR kd-Tree node to be modified.
    grid_id : `~yt.data_objects.grid_patch`
        A grid patch that contains the data spanned by this 
        kd-Tree leaf.
    leaf_l_corner: array_like, dimension 3
        The left corner of the volume spanned by this leaf.
    leaf_r_corner: array_like, dimension 3
        The right corner of the volume spanned by this leaf.
        
    Returns
    -------
    None
    """
    thisnode.grid = grid_id.id
    thisnode.l_corner = leaf_l_corner
    thisnode.r_corner = leaf_r_corner

class AMRKDTree(HomogenizedVolume):
    def __init__(self, pf,  l_max=None, le=None, re=None,
                 fields=None, no_ghost=False,
                 tree_type='domain',log_fields=None):
        r"""
        AMR kd-Tree object, a homogenized volume.

        Definition of the AMR kd-Tree object.  This is a method of
        volume homogenization that uses a modified kd-Tree structure
        to partition the AMR hierarchy.  The dividing nodes of the
        tree subdivide the volume into left and right children along a
        particular dimension.  The leaf nodes of the tree contain
        subvolumes that are covered by a single single grid at a
        single resolution, usually the maximum level for that volume
        unless `l_max` is otherwise specified.  The volume can then be
        traversed along an arbitrary direction based on comparisions
        with the dividing node position and split dimenstion.  

        Parameters
        ----------
        pf : `~yt.data_objects.StaticOutput`
            The parameter file to be kd-Tree partitioned.
        l_max : int, optional
            Maximum level to use in construction of kd-Tree. Default:
            None (all levels)
        le: array_like, optional
            Left edge to be be partitioned. Default: None (Domain Left
            Edge)
        re: array_like. optional
            Right edge to be partitioned.  Default: None (Domain Right
            Edge)
        fields: list of strings, optional
            Fields to be obtained when collecting leaf data.  Defualt:
            None (['Density']).
        log_fields: list of bool, optional
            Specifies which fields are to be taken the logarithm of
            before rendering.
        no_ghost: bool, optional
            Optimization option.  If True, homogenized bricks will
            extrapolate out from grid instead of interpolating from
            ghost zones that have to first be calculated.  This can
            lead to large speed improvements, but at a loss of
            accuracy/smoothness in resulting image.  The effects are
            less notable when the transfer function is smooth and
            broad. Default: False
        tree_type: string, optional
            Specifies the type of kd-Tree to be constructed/cast.
            There are three options, the default being 'domain'. Only
            affects parallel rendering.  'domain' is suggested.

            'domain' - Tree construction/casting is load balanced by
            splitting up the domain into the first N subtrees among N
            processors (N must be a power of 2).  Casting then
            proceeds with each processor rendering their subvolume,
            and final image is composited on the root processor.  The
            kd-Tree is never combined, reducing communication and
            memory overhead. The viewpoint can be changed without
            communication or re-partitioning of the data, making it
            ideal for rotations/spins.

            'breadth' - kd-Tree is first constructed as in 'domain',
            but then combined among all the subtrees.  Rendering is
            then split among N processors (again a power of 2), based
            on the N most expensive branches of the tree.  As in
            'domain', viewpoint can be changed without re-partitioning
            or communication.

            'depth' - kd-Tree is first constructed as in 'domain', but
            then combined among all subtrees.  Rendering is then load
            balanced in a back-to-front manner, splitting up the cost
            as evenly as possible.  If the viewpoint changes,
            additional data might have to be partitioned.  Is also
            prone to longer data IO times.  If all the data can fit in
            memory on each cpu, this can be the fastest option for
            multiple ray casts on the same dataset.


        Returns
        -------
        An AMR kd-Tree of the static output, of type `AMRKDTree`  

        Examples
        --------
        These are written in doctest format, and should illustrate how to
        use the function.  Use the variables 'pf' for the parameter file, 'pc' for
        a plot collection, 'c' for a center, and 'L' for a vector. 

        >>> from yt.utilities.amr_kdtree import AMRKDTree
        >>> volume = AMRKDTree(pf)
        yt         DEBUG      2010-11-08 21:35:40,873 Initializing data storage.
        yt         DEBUG      2010-11-08 21:35:40,873 Counting grids.
        yt         DEBUG      2010-11-08 21:35:40,874 Your data uses the annoying hardcoded path.
        yt         DEBUG      2010-11-08 21:35:40,876 Detected packed HDF5
        yt         DEBUG      2010-11-08 21:35:40,876 Setting up classes.
        yt         DEBUG      2010-11-08 21:35:40,877 Counting grids.
        yt         DEBUG      2010-11-08 21:35:40,877 Allocating arrays for 801 grids
        yt         DEBUG      2010-11-08 21:35:40,877 Parsing hierarchy.
        yt         INFO       2010-11-08 21:35:40,877 Getting the binary hierarchy
        yt         INFO       2010-11-08 21:35:40,885 Finished with binary hierarchy reading
        yt         DEBUG      2010-11-08 21:35:40,886 Constructing grid objects.
        yt         DEBUG      2010-11-08 21:35:40,903 Initializing data grid data IO
        yt         DEBUG      2010-11-08 21:35:40,904 Detecting fields.
        yt         DEBUG      2010-11-08 21:35:40,904 Adding unknown detected fields
        yt         DEBUG      2010-11-08 21:35:40,905 Setting up derived fields
        yt         DEBUG      2010-11-08 21:35:40,999 Re-examining hierarchy
        yt         INFO       2010-11-08 21:35:41,000 Making kd tree from le [ 0.  0.  0.] to [ 1.  1.  1.]
        yt         INFO       2010-11-08 21:35:42,451 Total of 5720 leafs
        yt         INFO       2010-11-08 21:35:42,519 [0000] Nodes 11439
        yt         INFO       2010-11-08 21:35:42,520 [0000] Cost is 314219
        yt         INFO       2010-11-08 21:35:42,520 [0000] Volume is 1.000000e+00
        >>> volume.volume
        1.0
        >>> volume.total_cost
        314219
        >>> volume.tree[0]
        {'cast_done': 0,
        'cost': 314219,
        'done': 0,
        'grid': None,
        'l_corner': array([ 0.,  0.,  0.]),
        'owner': 0,
        'r_corner': array([ 1.,  1.,  1.]),
        'split_ax': 1,
        'split_pos': 0.5}

        """

        self.pf = pf
        if nprocs > len(pf.h.grids):
            print('Parallel rendering requires that the number of \n \
            grids in the dataset is greater or equal to the number of \n \
            processors.  Reduce number of processors.')
            raise(KeyError)
        if fields is None: fields = ["Density"]
        self.no_ghost = no_ghost
        reduction_needed = {'domain':False,'depth':True,'breadth':True}
        self.tree_type = tree_type
        self.reduce_tree=reduction_needed[self.tree_type]
        self.bricks_loaded = False
        self.bricks = []

        self.fields = ensure_list(fields)
        if log_fields is not None:
            log_fields = ensure_list(log_fields)
        else:
            log_fields = [self.pf.field_info[field].take_log
                         for field in self.fields]
        self.log_fields = log_fields

        if l_max is None:
            self.l_max = self.pf.hierarchy.max_level+1
        else:
            self.l_max = na.min([l_max,self.pf.hierarchy.max_level+1])

        if le is None:
            self.domain_left_edge = pf.domain_left_edge
        else:
            self.domain_left_edge = na.clip(na.array(le),pf.domain_left_edge, pf.domain_right_edge)
        if re is None:
            self.domain_right_edge = pf.domain_right_edge
        else:
            self.domain_right_edge = na.clip(na.array(re),pf.domain_left_edge, pf.domain_right_edge)

        self.my_l_corner = self.domain_left_edge
        self.my_r_corner = self.domain_right_edge

        mylog.info('Making kd tree from le %s to %s'% (self.domain_left_edge, self.domain_right_edge))
        root_grids = pf.hierarchy.get_levels().next()

        root_l_data = na.array([grid.LeftEdge for grid in root_grids])
        root_r_data = na.array([grid.RightEdge for grid in root_grids])
        root_we_want = na.all(root_l_data < self.my_r_corner,axis=1)*\
                       na.all(root_r_data > self.my_l_corner,axis=1)
        
        root_grids = root_grids[root_we_want]

        # Build the kd-Tree
        nodes = self.__build(root_grids, None, self.domain_left_edge, self.domain_right_edge)

        # Set up kd-Tree Dictionary
        mytree = {}
        if self.reduce_tree:
            for node in nodes:
                if node.owner is my_rank:
                    mytree[node.id] = {'l_corner':node.l_corner, 'r_corner':node.r_corner,
                                       'grid':node.grid, 'split_ax':node.split_ax, 'split_pos':node.split_pos, 'owner':node.owner}
        else:
            for node in nodes:
                mytree[node.id] = {'l_corner':node.l_corner, 'r_corner':node.r_corner,
                                   'grid':node.grid, 'split_ax':node.split_ax, 'split_pos':node.split_pos, 'owner':node.owner}

        # Merge the kd-trees (only applies in parallel)
        self.merge_trees(mytree)
        # Add properties to leafs/nodes
        self.rebuild_tree_links()
        # Calculate the total cost of the render
        self.total_cost = self.count_cost()
        # Calculate the total volume spanned by the tree
        self.volume = self.count_volume()
        mylog.info('[%04i] Nodes %d' % (my_rank,len(self.tree)))
        mylog.info('[%04i] Cost is %d' % (my_rank,self.total_cost))
        mylog.info('[%04i] Volume is %e' % (my_rank,self.volume)) 
        
    def merge_trees(self, mytree):
        if nprocs > 1 and self.reduce_tree:
            self.tree = self._mpi_joindict(mytree)
        else:
            self.tree = mytree

    def get_bricks(self):
        r"""Preload the bricks into the kd-Tree

        Traverses the tree, gets the vertex centered data, and
        attaches partitioned grids to the kd-Tree structure.
        
        Parameters
        ----------
        None

        Returns
        ----------
        None
        
        """
        if self.bricks_loaded: return
        current_saved_grids = []
        current_vcds = []
        
        for current_node in self.tree.itervalues():
            if current_node['grid'] is None: continue

            if current_node['grid'] in current_saved_grids:
                dds = current_vcds[current_saved_grids.index(current_node['grid'])]
            else:
                dds = []
                for i,field in enumerate(self.fields):
                    vcd = current_node['grid'].get_vertex_centered_data(field,smoothed=True,no_ghost=self.no_ghost).astype('float64')
                    if self.log_fields[i]: vcd = na.log10(vcd)
                    dds.append(vcd)
                current_saved_grids.append(current_node['grid'])
                current_vcds.append(dds)

            data = [d[current_node['li'][0]:current_node['ri'][0]+1,
                      current_node['li'][1]:current_node['ri'][1]+1,
                      current_node['li'][2]:current_node['ri'][2]+1].copy() for d in dds]
            
            current_node['brick'] = PartitionedGrid(current_node['grid'].id, len(self.fields), data,
                                                    current_node['l_corner'].copy(), 
                                                    current_node['r_corner'].copy(), 
                                                    current_node['dims'].astype('int64'))
            self.bricks.append(current_node['brick'])
        del current_saved_grids, current_vcds
        self.bricks_loaded = True
        
    def rebuild_tree_links(self):
        r"""Sets properties of the kd-Tree dictionary

        For each node, sets the `cost`, `done`, `cast_done`, and `owner` keys.
        
        Parameters
        ----------
        None

        Returns
        ----------
        None
        
        """
        num_leafs = 0
        for this_id, node in self.tree.iteritems():
            node['cost']=0
            if node['grid'] is not None:
                # I'm a leaf!
                self.set_leaf_props(node)
                num_leafs+=1
            node['done'] = 0
            node['cast_done'] = 0
            if this_id < nprocs-1:
                node['owner'] = (this_id+1 - 2**(len(na.binary_repr(this_id+1))-1))*\
                                2**(len(na.binary_repr(nprocs))-len(na.binary_repr(this_id+1)))
        mylog.info('Total of %i leafs' % num_leafs)

    def set_leaf_props(self,thisnode):
        r"""Given a leaf, gathers grid, indices, dimensions, and cost properties.

        Parameters
        ----------
        None

        Returns
        ----------
        None
        
        """
        thisnode['grid'] = self.pf.hierarchy.grids[thisnode['grid'] - 1]

        dds = thisnode['grid'].dds
        gle = thisnode['grid'].LeftEdge
        gre = thisnode['grid'].RightEdge
        thisnode['li'] = ((thisnode['l_corner']-gle)/dds).astype('int32')
        thisnode['ri'] = ((thisnode['r_corner']-gle)/dds).astype('int32')
        thisnode['dims'] = (thisnode['ri'] - thisnode['li']).astype('int32')
        # Here the cost is actually inversely proportional to 4**Level (empirical)
        thisnode['cost'] = (na.prod(thisnode['dims'])/4.**thisnode['grid'].Level).astype('int64')
        # Here is the old way
        # thisnode.cost = na.prod(thisnode.dims).astype('int64')
        

    def count_cost(self):
        r"""Counts the cost of the entire tree, while filling in branch costs.

        Parameters
        ----------
        None

        Returns
        ----------
        Total cost of rendering the kd-Tree

        At completion, each node in the kd-Tree carries the total cost
        of all branches and leaves it contains.
        
        """
        ids = self.tree.keys()
        ids.sort()
        total_cost = 0
        for this_id in ids[::-1]:
            if self.tree[this_id]['grid'] is not None:
                total_cost += self.tree[this_id]['cost']
            try:
                self.tree[_parent_id(this_id)]['cost'] += self.tree[this_id]['cost']
            except:
                pass
        return total_cost
        
    def count_volume(self):
        r"""Calculates the volume of the kd-Tree

        Parameters
        ----------
        None

        Returns
        ----------
        Total volume of the tree.
        
        """
        v = 0.0
        for node in self.tree.itervalues():
            if node['grid'] is not None:
                v += na.prod(node['r_corner'] - node['l_corner'])
        return v

    def reset_cast(self):
        r"""Resets `cast_done` to 0"""
        for node in self.tree.itervalues():
            node['cast_done'] = 0

    def __build(self, grids, parent, l_corner, r_corner):
        r"""Builds the AMR kd-Tree

        Parameters
        ----------
        grids: array_like
            Array of grids that cover the volume to be decomposed into
            the kd-Tree
        parent: ~yt.data_objects.grid_patch
            The parent grid that covers the volume.  Can be None if
            the volume is not contained by a single grid.
        l_corner: array_like
            The left corner of the volume to be decomposed.
        r_corner: array_like
            The right corner of the volume to be decomposed
            
        Returns
        ----------
        An array of kd-Tree nodes that make up the AMR kd-Tree
        
        """
        nodes = [MasterNode()]
        
        head_node = nodes[0]
        current_node = nodes[0]
        current_node.grids = grids
        current_node.l_corner = l_corner
        current_node.r_corner = r_corner
        current_node.owner = my_rank
        current_node.id = 0
        par_tree_depth = long(na.log2(nprocs))

        while(not head_node.done):

            if current_node.done:
                current_node = current_node.parent
                # print 'I am going to work on my parent because I am done'
                continue

            if current_node.left_children is not None and current_node.right_children is not None:
                if current_node.left_children.done and current_node.right_children.done:
                    current_node.done = 1
                    # print 'I am done with this node because my children are done.'
                    continue
                elif not current_node.left_children.done:
                    if ((current_node.left_children.id + 1)>>par_tree_depth) == 1:
                        # There are nprocs nodes that meet this criteria
                        if (current_node.left_children.id+1-nprocs) is my_rank:
                            # I own this shared node
                            current_node = current_node.left_children
                            current_node.owner = my_rank
                            
                            # print '[%04i]'%my_rank,'I am building the kD-tree in the region:',\
                            # current_node.l_corner,\
                            # current_node.r_corner,\
                            # na.prod(current_node.r_corner-current_node.l_corner),\
                            # len(current_node.grids), current_node.id
                            self.my_l_corner = current_node.l_corner
                            self.my_r_corner = current_node.r_corner
                            continue
                        else:
                            # I don't own this shared node
                            current_node.left_children.done = 1
                            current_node.left_children.owner = my_rank-1
                            continue
                    else:
                        # Otherwise move to the left child
                        current_node = current_node.left_children
                        current_node.owner = my_rank
                        continue
                elif not current_node.right_children.done:
                    if ((current_node.right_children.id + 1)>>par_tree_depth) == 1:
                        # There are nprocs nodes that meet this criteria
                        if (current_node.right_children.id+1-nprocs) is my_rank:
                            # I own this shared node
                            current_node = current_node.right_children
                            current_node.owner = my_rank

                            # print '[%04i]'%my_rank,'I am building the kD-tree in the region:',\
                            # current_node.l_corner,\
                            # current_node.r_corner,\
                            # na.prod(current_node.r_corner-current_node.l_corner),\
                            # len(current_node.grids), current_node.id
                            self.my_l_corner = current_node.l_corner
                            self.my_r_corner = current_node.r_corner

                            continue
                        else:
                            # I don't own this shared node
                            current_node.right_children.done = 1
                            current_node.right_children.owner = my_rank+1
                            continue
                    else:
                        # Otherwise move to the right child
                        current_node = current_node.right_children
                        current_node.owner = my_rank
                        continue
                
            
            # If we are in a single grid
            if len(current_node.grids) is 1:
                thisgrid = current_node.grids[0]
                # If we are in the specified domain
                if (thisgrid.LeftEdge[0] <= current_node.l_corner[0]) and (thisgrid.RightEdge[0] >= current_node.r_corner[0]) and \
                   (thisgrid.LeftEdge[1] <= current_node.l_corner[1]) and (thisgrid.RightEdge[1] >= current_node.r_corner[1]) and \
                   (thisgrid.LeftEdge[2] <= current_node.l_corner[2]) and (thisgrid.RightEdge[2] >= current_node.r_corner[2]):
                    children = []
                    # Check if we have children and have not exceeded l_max
                    if len(thisgrid.Children) > 0 and thisgrid.Level < self.l_max:
                        children = self.pf.hierarchy.grids[na.array([child.id - 1 for child in thisgrid.Children])]

                        child_l_data = na.array([child.LeftEdge for child in children])
                        child_r_data = na.array([child.RightEdge for child in children])
                        children_we_want = na.all(child_l_data < current_node.r_corner,axis=1)*\
                                           na.all(child_r_data > current_node.l_corner,axis=1)

                        children = children[children_we_want]
                        del child_l_data, child_r_data, children_we_want
                    # Continue with either our children or make a leaf out of ourself.
                    if len(children) is 0:
                        set_leaf(current_node, thisgrid, current_node.l_corner, current_node.r_corner)
                        current_node.done = 1
                        # print 'My single grid covers the rest of the volume, and I have no children'
                        del children
                        continue
                    else:
                        current_node.grids = children
                        current_node.parent_grid = thisgrid
                        # print 'My single grid covers the rest of the volume, and I children, about to iterate on them'
                        del children
                        continue
            # If we don't have any grids, this volume belongs to the parent        
            if len(current_node.grids) is 0:
                set_leaf(current_node, current_node.parent_grid, current_node.l_corner, current_node.r_corner)
                current_node.done = 1
                # print 'This volume does not have a child grid, so it belongs to my parent!'
                continue

            # Get the left edges for each child
            l_data = na.array([child.LeftEdge for child in current_node.grids])
            r_data = na.array([child.RightEdge for child in current_node.grids])

            # Split along the best dimension
            best_dim = 0
            best_choices = []

            # This helps a lot in terms of speed, but in some crazy situation could cause problems.
            if l_data.shape[0] > 20:
                best_dim = na.argmax(current_node.r_corner - current_node.l_corner)
                choices = na.concatenate(([current_node.l_corner[best_dim]],l_data[:,best_dim],r_data[:,best_dim],
                                          [current_node.r_corner[best_dim]]))
                choices = choices[choices >= current_node.l_corner[best_dim]]
                best_choices = choices[choices <= current_node.r_corner[best_dim]]
                best_choices.sort()
            else:
                for d in range(3):
                    choices = na.unique(na.clip(na.concatenate(([current_node.l_corner[d]],l_data[:,d],r_data[:,d],
                                                                [current_node.r_corner[d]])),
                                                current_node.l_corner[d],current_node.r_corner[d]))
                    if len(choices) > len(best_choices):
                        best_choices = choices
                        best_dim = d

            split = best_choices[len(best_choices)/2]

            less_ids = na.nonzero(l_data[:,best_dim] < split)[0]
            greater_ids = na.nonzero(split < r_data[:,best_dim])[0]

            current_node.split_ax = best_dim
            current_node.split_pos = split

            # print 'child has %i grids' % sum(less_ids)
            nodes.append(MasterNode())
            current_node.left_children = nodes[-1]
            current_node.left_children.id = _lchild_id(current_node.id)
            current_node.left_children.parent = current_node
            current_node.left_children.parent_grid = current_node.parent_grid
            current_node.left_children.grids = current_node.grids[less_ids]
            current_node.left_children.l_corner = current_node.l_corner
            current_node.left_children.r_corner = corner_bounds(best_dim, split, current_right=current_node.r_corner)

            # print 'child has %i grids' % sum(greater_ids)
            nodes.append(MasterNode())
            current_node.right_children = nodes[-1]
            current_node.right_children.id = _rchild_id(current_node.id)
            current_node.right_children.parent = current_node
            current_node.right_children.parent_grid = current_node.parent_grid
            current_node.right_children.grids = current_node.grids[greater_ids]
            current_node.right_children.l_corner = corner_bounds(best_dim, split, current_left=current_node.l_corner)
            current_node.right_children.r_corner = current_node.r_corner

            del l_data, r_data, best_dim, best_choices, split, less_ids, greater_ids, choices

        return nodes


    def traverse(self, back_center, front_center, start_id):
        r"""Traverses the kd-Tree, returning a list of partitioned grids.

        Parameters
        ----------
        back_center: array_like
            Position of the back center from which to start moving forward.
        front_center: array_like
            Position of the front center to which the traversal progresses.
        start_id: int
            First kd-Tree node to begin with

        Returns
        ----------
        An array of partitioned grids, ordered from the back_center to
        the front_center.
        
        """

        viewpoint = front_center - back_center
        current_id = start_id
        tree = self.tree
        head_node = tree[current_id]
        current_node = tree[current_id]
        if self.tree_type is 'depth':
            cast_all = False
        else:
            cast_all = True
        total_cells = 0
        my_total = 0
        p_grids = []

        self.current_saved_grids = []
        self.current_vcds = []

        while(True):

            current_node = tree[current_id]
            
            if head_node['cast_done'] is 1:
                break

            if current_node['cast_done'] is 1:
                current_id = _parent_id(current_id)
                # print 'I am going to work on my parent because I am done'
                continue

            if current_node['grid'] is None:
                if _lchild_id(current_id) not in tree or tree[_lchild_id(current_id)]['cast_done'] == 1:
                    if _rchild_id(current_id) not in tree or tree[_rchild_id(current_id)]['cast_done'] == 1:
                        current_node['cast_done'] = 1
                        # print 'I am done with this node because my children are done.'
                        continue

            if current_node['grid'] is None:
                # print current_id, current_node, viewpoint, current_node['split_ax'], current_node['split_pos']
                if viewpoint[current_node['split_ax']] <= current_node['split_pos']:
                    if _lchild_id(current_id) in tree:
                        if tree[_lchild_id(current_id)]['cast_done'] == 0:
                            current_id = _lchild_id(current_id)
                            # print 'I am going to work on my left child'
                            continue
                    if _rchild_id(current_id) in tree:
                        if tree[_rchild_id(current_id)]['cast_done'] == 0:
                            current_id = _rchild_id(current_id)
                            # print 'I am going to work on my right child'
                            continue
                else:
                    if _rchild_id(current_id) in tree:
                        if tree[_rchild_id(current_id)]['cast_done'] == 0:
                            current_id = _rchild_id(current_id)
                            # print 'I am going to work on my right child'
                            continue
                    if _lchild_id(current_id) in tree:
                        if tree[_lchild_id(current_id)]['cast_done'] == 0:
                            current_id = _lchild_id(current_id)
                            # print 'I am going to work on my left child'
                            continue
            else:
                if (((total_cells >= 1.0*my_rank*self.total_cost/nprocs) and
                     (total_cells < 1.0*(my_rank+1)*self.total_cost/nprocs)) or cast_all):
                    # print ' I am actually casting'
                    if current_node.get('brick') is None:
                        #current_node_time = time.time()
                        if current_node['grid'] in self.current_saved_grids:
                            dds = self.current_vcds[self.current_saved_grids.index(current_node['grid'])]
                        else:
                            dds = []
                            for i,field in enumerate(self.fields):
                                vcd = current_node['grid'].get_vertex_centered_data(field,no_ghost=self.no_ghost).astype('float64')
                                if self.log_fields[i]: vcd = na.log10(vcd)
                                dds.append(vcd)
                            self.current_saved_grids.append(current_node['grid'])
                            self.current_vcds.append(dds)

                        data = [d[current_node['li'][0]:current_node['ri'][0]+1,
                                  current_node['li'][1]:current_node['ri'][1]+1,
                                  current_node['li'][2]:current_node['ri'][2]+1].copy() for d in dds]
                        
                        current_node['brick'] = PartitionedGrid(current_node['grid'].id, len(self.fields), data,
                                                                current_node['l_corner'].copy(), 
                                                                current_node['r_corner'].copy(), 
                                                                current_node['dims'].astype('int64'))
                    p_grids.append(current_node)
                    my_total += current_node['cost']
                total_cells += current_node['cost']
                current_node['cast_done'] = 1
        return p_grids

    def initialize_source(self):
        if self.tree_type is 'domain':
            self.get_bricks()
        
    def kd_ray_cast(self,image, tfp, vector_plane, back_center, front_center):
        r"""Traverses the kd-Tree, casting the partitioned grids from back to
            front.

        Given an image, transfer function, vector plane, back and
        front centers, ray-cast using the kd-Tree structure.

        Parameters
        ----------
        image: na.array
            Image plane to contain resulting ray cast.
        tfp: ~yt.utilities.amr_utils.TransferFunctionProxy
            Transfer function to be used in rendering.
        vector_plane: ~yt.utilities.amr_utils.VectorPlane
            Vector plane object to be used during ray casting
        back_center: array_like
            Position of the back center from which to start moving forward.
        front_center: array_like
            Position of the front center to which the traversal progresses.
            
        Returns
        ----------
        image: na.array
            An rgb array of the resulting rendering.
        
        See Also
        ----------
        yt.visualization.volume_rendering.camera
        
        """
        if self.tree is None: 
            print 'No KD Tree Exists'
            return
        self.total_cells = 0
        self.image = image
        if self.tree_type is 'domain':
            depth = int(na.log2(nprocs))
            start_id = 0

            rt1 = time()
            mylog.info('About to cast')

            pbar = get_pbar("Ray casting",self.total_cost)
            total_cells = 0
            for brick in self.traverse(back_center, front_center, start_id):
                brick['brick'].cast_plane(tfp, vector_plane)
                total_cells += brick['cost']
                pbar.update(total_cells)
            pbar.finish()
        
            mylog.info('I am done with my rendering after %e seconds' % (time()-rt1)) 
            self.reduce_tree_images(self.tree, front_center)
            mylog.info('Done in kd_ray_cast') 

        elif self.tree_type is 'breadth':
            if self.reduce_tree is False and nprocs > 1:
                mylog.error("Breadth-first rendering requires keyword volume_rendering(merge_kd_tree=True) ")
                mylog.error("Perhaps try the cast_type 'domain'")
                raise(KeyError)
                    
            self.clean_owners(self.tree)
            rt1 = time()
            ids = [0]
            costs = [self.tree[0]['cost']]

            new_ids, new_costs = self.split_tree(ids,costs,nprocs)
            for i,this_id in enumerate(new_ids):
                self.tree[this_id]['owner'] = i
                if my_rank == self.tree[this_id]['owner']:
                    my_node = this_id

                    pbar = get_pbar("Ray casting",self.total_cost)
                    total_cells = 0
                    for brick in self.traverse(back_center, front_center, this_id):
                        brick['brick'].cast_plane(tfp, vector_plane)
                        total_cells += na.prod(brick['cost'])
                        pbar.update(total_cells)
                    pbar.finish()
        
            print '[%04d] I am done with my rendering after %e seconds' % (my_rank, time()-rt1) 
            self.breadth_reduce_tree(0,front_center-back_center)
            final_owner = self.tree[0]['owner']
            if final_owner != 0:
                if final_owner == my_rank:
                    PT._send_array(self.image.ravel(), 0, tag=final_owner)
                if my_rank == 0:
                    self.image = PT._recv_array(final_owner, tag=final_owner).reshape(
                        (self.image.shape[0],self.image.shape[1],self.image.shape[2]))
            print 'Done in kd_ray_cast' 

        elif self.tree_type is 'depth':
            if self.reduce_tree is False and nprocs > 1:
                mylog.error("Depth-first rendering requires keyword volume_rendering(merge_kd_tree=True) ")
                mylog.error("Perhaps try the cast_type 'domain'")
                raise(KeyError)

            mylog.info('About to cast')
            rt1 = time()
            pbar = get_pbar("Ray casting",self.total_cost)
            total_cells = 0
            for brick in self.traverse(back_center, front_center, 0):
                brick['brick'].cast_plane(tfp, vector_plane)
                total_cells += na.prod(brick['cost'])
                pbar.update(total_cells)
            pbar.finish()
                
            mylog.info('I am done with my rendering after %e seconds' % (time()-rt1)) 
            im = self._binary_tree_reduce(self.image)
            self.image = im
            mylog.info('Done in kd_ray_cast')
            self._barrier()
        return self.image

    def reduce_tree_images(self, tree, viewpoint):
        rounds = int(na.log2(nprocs))
        my_node = 2**rounds - 1 + my_rank
        for thisround in range(rounds,0,-1):
            #print my_rank, 'my node', my_node
            parent_id = _parent_id(my_node)
            parent = tree[parent_id]
            #print parent['split_ax'], parent['split_pos']
            if viewpoint[parent['split_ax']] <= parent['split_pos']:
                front = tree[_lchild_id(parent_id)]
                back = tree[_rchild_id(parent_id)]
            else:
                front = tree[_rchild_id(parent_id)]
                back = tree[_lchild_id(parent_id)]

            # Send the images around
            if front['owner'] == my_rank:
                if front['owner'] == parent['owner']:
                    #print my_rank, 'receiving image from ',back['owner']
                    arr2 = PT._recv_array(back['owner'], tag=back['owner']).reshape(
                        (self.image.shape[0],self.image.shape[1],self.image.shape[2]))
                    for i in range(3):
                        # This is the new way: alpha corresponds to opacity of a given
                        # slice.  Previously it was ill-defined, but represented some
                        # measure of emissivity.
                        #print arr2.shape
                        #                ta = (1.0 - arr2[:,:,i+3])
                        ta = (1.0 - na.sum(self.image,axis=2))
                        ta[ta<0.0] = 0.0 
                        self.image[:,:,i  ] = self.image[:,:,i  ] + ta * arr2[:,:,i  ]
                else:
                    mylog.info('Reducing image.  You have %i rounds to go in this binary tree' % thisround)
                    #print my_rank, 'sending my image to ',back['owner']
                    PT._send_array(self.image.ravel(), back['owner'], tag=my_rank)

                
            if back['owner'] == my_rank:
                if front['owner'] == parent['owner']:
                    #print my_rank, 'sending my image to ',front['owner']
                    PT._send_array(self.image.ravel(), front['owner'], tag=my_rank)
                else:
                    mylog.info('Reducing image.  You have %i rounds to go in this binary tree' % thisround)
                    #print my_rank, 'receiving image from ',front['owner']
                    arr2 = PT._recv_array(front['owner'], tag=front['owner']).reshape(
                        (self.image.shape[0],self.image.shape[1],self.image.shape[2]))
                    for i in range(3):
                        # This is the new way: alpha corresponds to opacity of a given
                        # slice.  Previously it was ill-defined, but represented some
                        # measure of emissivity.
                        # print arr2.shape
                        # ta = (1.0 - arr2[:,:,i+3])
                        ta = (1.0 - na.sum(arr2,axis=2))
                        ta[ta<0.0] = 0.0 
                        self.image[:,:,i  ] = arr2[:,:,i  ] + ta * self.image[:,:,i  ]
                        # image[:,:,i+3] = arr2[:,:,i+3] + ta * image[:,:,i+3]
            # Set parent owner to back owner
            my_node = (my_node-1)>>1

    def breadth_reduce_tree(self,this_id, viewpoint):
        my_rank = self._mpi_get_rank()
        tree = self.tree
        if tree[this_id]['owner'] == my_rank:
            return
        if viewpoint[tree[this_id]['split_ax']] <= tree[this_id]['split_pos']:
            # Stuff on the right is further away than left, reduce left on top of right.
            front = _lchild_id(this_id)
            back = _rchild_id(this_id)
        else:
            # Stuff on the right is closer than left, reduce right on top of left.
            front = _rchild_id(this_id)
            back = _lchild_id(this_id)

    
        # reduce the children first
        if tree[back]['owner'] == -1:
            print 'back owner is -1, reducing it'
            self.breadth_reduce_tree(back, viewpoint)
        if tree[front]['owner'] == -1:
            print 'front owner is -1, reducing it'
            self.breadth_reduce_tree(front,viewpoint)

        # Send the images around
        if tree[front]['owner'] == my_rank:
            print my_rank, 'sending my image to ',tree[back]['owner']
            PT._send_array(self.image.ravel(), tree[back]['owner'], tag=my_rank)
        if tree[back]['owner'] == my_rank:
            print my_rank, 'receiving image from ',tree[front]['owner']
            arr2 = PT._recv_array(tree[front]['owner'], tag=tree[front]['owner']).reshape(
                (self.image.shape[0],self.image.shape[1],self.image.shape[2]))
            for i in range(3):
                # print arr2.shape
                # ta = (1.0 - arr2[:,:,i+3])
                ta = (1.0 - na.sum(arr2,axis=2))
                ta[ta<0.0] = 0.0 
                self.image[:,:,i  ] = arr2[:,:,i  ] + ta * self.image[:,:,i  ]
                # image[:,:,i+3] = arr2[:,:,i+3] + ta * image[:,:,i+3]
        # Set parent owner to back owner
        tree[this_id]['owner'] = tree[back]['owner']
        return

    def _binary_tree_reduce(self, arr, root = 0):
        self._barrier()
        myrank = self._mpi_get_rank()
        nprocs = self._mpi_get_size()
        depth = 0
        procs_with_images = range(0,nprocs,2**depth)
        
        while(len(procs_with_images) > 1):
            if not myrank in procs_with_images:
                break
            # Even - receive then composite
            if procs_with_images.index(myrank)%2 == 0:
                arr2 = PT._recv_array(myrank+2**depth, tag=myrank+2**depth).reshape(
                    (arr.shape[0],arr.shape[1],arr.shape[2]))

                for i in range(3):
                    # This is the new way: alpha corresponds to opacity of a given
                    # slice.  Previously it was ill-defined, but represented some
                    # measure of emissivity.
                    ta = (1.0 - na.sum(arr2,axis=2))
                    ta[ta<0.0] = 0.0 
                    arr[:,:,i  ] = arr2[:,:,i  ] + ta * arr[:,:,i  ]
                    del ta
            # Odd - send
            else:
                PT._send_array(arr.ravel(), myrank-2**depth, tag=myrank)
            depth += 1
            procs_with_images = range(0,nprocs,2**depth) 
        self._barrier()
        return arr

    def split_tree(self,ids, costs, n_trees):
        tree = self.tree
        if len(ids) == n_trees:
            return ids, costs

        expensive_index = na.argmax(costs)
        expensive_tree = ids[expensive_index]
        
        if tree[ids[expensive_index]]['grid'] is None:
            tree[ids[expensive_index]]['owner'] = -1
            to_split = ids.pop(expensive_index)
            del costs[expensive_index]
            lchild_id = _lchild_id(expensive_tree)
            ids.append(lchild_id)
            costs.append(tree[lchild_id]['cost'])
            
            rchild_id = _rchild_id(expensive_tree)
            ids.append(rchild_id)
            costs.append(tree[rchild_id]['cost'])
            new_ids, new_costs = self.split_tree(ids, costs, n_trees)
            
        else:
            expensive_leaf = ids.pop(expensive_index)
            leaf_cost = costs.pop(expensive_index)
            new_ids, new_costs = self.split_tree(ids, costs, n_trees-1)
            new_ids.insert(expensive_tree,expensive_leaf)
            new_costs.insert(expensive_tree,leaf_cost)

        return new_ids, new_costs

    def clean_owners(self,nodes):
        r"""Sets all kd-Tree node `owner` to -1"""
        for node in nodes.itervalues():
            node['owner'] = -1
 
    def store_kd_bricks(self, fn=None):
        if fn is None:
            fn = '%s_kd_bricks.h5'%self.pf
        if my_rank != 0:
            PT._recv_array(my_rank-1, tag=my_rank-1)
        f = h5py.File(fn,"a")
        for i, node in self.tree.iteritems():
            if 'brick' in node:
                for fi,field in enumerate(self.fields):
                    try:
                        f.create_dataset("/brick_%s_%s" % (hex(i),field),
                                         data = node['brick'].my_data[fi].astype('float64'))
                    except:
                        pass
        f.close()
        if my_rank != (nprocs-1):
            PT._send_array([0],my_rank+1, tag=my_rank)
        
    def load_kd_bricks(self,fn=None):
        if fn is None:
            fn = '%s_kd_bricks.h5' % self.pf
        if my_rank != 0:
            PT._recv_array(my_rank-1, tag=my_rank-1)
        try:
            f = h5py.File(fn,"r")
            for i, node in self.tree.iteritems():
                if node['grid'] is not None:
                        data = [f["brick_%s_%s" %
                                           (hex(i), field)][:].astype('float64') for field in self.fields]
                        node['brick'] = PartitionedGrid(node['grid'].id, len(self.fields), data,
                                                        node['l_corner'].copy(), 
                                                        node['r_corner'].copy(), 
                                                        node['dims'].astype('int64'))
            self.bricks_loaded=True
            f.close()
        except:
            pass
        if my_rank != (nprocs-1):
            PT._send_array([0],my_rank+1, tag=my_rank)

    def load_tree(self,fn):
        f = h5py.File(fn,"r")
        kd_ids = f["/kd_ids"][:]
        kd_l_corners = f['/left_edges'][:]
        kd_r_corners = f['/right_edges'][:]
        kd_grids = f['/grids'][:]
        kd_split_axs = f['/split_axs'][:]
        kd_split_pos = f['/split_pos'][:]
        kd_owners = f['/kd_owners'][:]

        mytree = {}
        for i,this_id in enumerate(kd_ids):
            mytree[this_id] = {'l_corner':kd_l_corners[i], 'r_corner':kd_r_corners[i],
                               'split_ax':kd_split_axs[i], 'split_pos':kd_split_pos[i], 'owner':kd_owners[i]}
            if kd_grids[i] == -1:
                mytree[this_id]['grid'] = None
                mytree[this_id]['brick'] = None
            else:
                mytree[this_id]['grid'] = kd_grids[i]
                mytree[this_id]['brick'] = [f["/bricks/brick_%s/%s" % (i, field)][:] for field in self.fields]
                mytree[this_id]['split_ax'] = None
                mytree[this_id]['split_pos'] = None
        f.close()

        self.tree = mytree


    def store_tree(self,fn):
        f = h5py.File(fn,"w")
        Nkd = len(self.tree)
        kd_l_corners = na.zeros( (Nkd, 3), dtype='float64')
        kd_r_corners = na.zeros( (Nkd, 3), dtype='float64')
        kd_grids = na.zeros( (Nkd) )
        kd_split_axs = na.zeros( (Nkd), dtype='int32')
        kd_split_pos = na.zeros( (Nkd), dtype='float64')
        kd_owners = na.zeros( (Nkd), dtype='int32')
        f.create_group("/bricks")
        for i, tree_item in enumerate(self.tree.iteritems()):
            kdid = tree_item[0]
            node = tree_item[1]
            kd_l_corners[i,:] = node['l_corner']
            kd_r_corners[i,:] = node['r_corner']
            if node['grid'] is None:
                kd_grids[i] = -1
                kd_split_axs[i] = node['split_ax']
                kd_split_pos[i] = node['split_pos']
            else:
                kd_grids[i] = node['grid'].id
                kd_split_axs[i] = -1
                kd_split_pos[i] = 0.0
            
            kd_owners[i] = node['owner']
            if 'brick' in node:
                f.create_group("/bricks/brick_%08i" % i)
                for fi,field in enumerate(self.fields):
                    f.create_dataset("/bricks/brick_%08i/%s" % (i,field),
                                     data = node['brick'].my_data[fi])
        f.create_dataset("/left_edges",data=kd_l_corners)
        f.create_dataset("/right_edges",data=kd_r_corners)
        f.create_dataset("/grids",data=kd_grids)
        f.create_dataset("/split_axs",data=kd_split_axs)
        f.create_dataset("/split_pos",data=kd_split_pos)
        f.create_dataset("/kd_owners",data=kd_owners)
        f.close()
