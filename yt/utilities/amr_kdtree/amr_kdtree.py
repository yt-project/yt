"""
AMR kD-Tree Framework

Authors: Samuel Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
Wil St. Charles <fallen751@gmail.com>
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
from yt.funcs import *
from yt.visualization.volume_rendering.grid_partitioner import HomogenizedVolume
from yt.visualization.image_writer import write_image, write_bitmap
from yt.utilities.amr_utils import kdtree_get_choices
from yt.utilities._amr_utils.grid_traversal import PartitionedGrid
from yt.utilities.performance_counters import yt_counters, time_function
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import ParallelAnalysisInterface 
import matplotlib.pylab as pl
from copy import deepcopy
from yt.config import ytcfg
from time import time
import h5py

def corner_bounds(split_dim, split, current_left = None, current_right = None):
    r"""
    Given a kd-Tree split dimension and position and bound to be
    modified, returns the new bound.

    A simple function that replaces the `split_dim` dimension of the
    current left or right bound with `split`. Left or Right bound is
    chosen by specifying the `current_left` or `current_right`.
    """
    if(current_left is not None):
        new_left = current_left.copy()
        new_left[split_dim] = split
        return new_left
    elif(current_right is not None):
        new_right = current_right.copy()
        new_right[split_dim] = split
        return new_right

def _lchild_id(id): return (id<<1) + 1
def _rchild_id(id): return (id<<1) + 2
def _parent_id(id): return (id-1)>>1

steps = na.array([[-1, -1, -1],
                  [-1, -1,  0],
                  [-1, -1,  1],
                  [-1,  0, -1],
                  [-1,  0,  0],
                  [-1,  0,  1],
                  [-1,  1, -1],
                  [-1,  1,  0],
                  [-1,  1,  1],
                  
                  [ 0, -1, -1],
                  [ 0, -1,  0],
                  [ 0, -1,  1],
                  [ 0,  0, -1],
                  # [ 0,  0,  0],
                  [ 0,  0,  1],
                  [ 0,  1, -1],
                  [ 0,  1,  0],
                  [ 0,  1,  1],
                  
                  [ 1, -1, -1],
                  [ 1, -1,  0],
                  [ 1, -1,  1],
                  [ 1,  0, -1],
                  [ 1,  0,  0],
                  [ 1,  0,  1],
                  [ 1,  1, -1],
                  [ 1,  1,  0],
                  [ 1,  1,  1]
                  ])


class MasterNode(object):
    r"""
    A MasterNode object is the building block of the AMR kd-Tree.
    Used during the construction to act as both dividing nodes and
    leaf nodes.
    """
    def __init__(self,my_id=None,
                 parent=None,
                 parent_grid=None,
                 grids=None,
                 l_corner=None,
                 r_corner=None):
        
        self.grids = grids
        self.parent = parent
        self.parent_grid = parent_grid
        self.l_corner = l_corner
        self.r_corner = r_corner

        self.split_ax = None
        self.split_pos = None

        self.left_child = None
        self.right_child = None
        self.cost = 0
        # self.owner = -1

        self.id = my_id
        
        self.grid = None
        self.brick = None
        self.li = None
        self.ri = None
        self.dims = None

        # self.done = 0
        # self.cast_done = 0

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
    del thisnode.grids, thisnode.parent_grid, thisnode.split_ax, thisnode.split_pos

class AMRKDTree(HomogenizedVolume):
    def __init__(self, pf,  l_max=None, le=None, re=None,
                 fields=None, no_ghost=False,
                 tree_type='domain',log_fields=None, merge_trees=False):
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
        merge_trees: bool, optional
            If True, the distributed kD-tree can be merged
            together in an allgather-like fashion.  This should not be
            used for parallel rendering, as it will cause all
            processors to render all bricks.  This is primarily useful
            for applications that are not domain decomposed but still
            want to build the kD-tree in parallel. Default:False
        

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
        ParallelAnalysisInterface.__init__(self)
        self.current_split_dim = 0

        self.pf = pf
        self.sdx = self.pf.h.get_smallest_dx()
        self._id_offset = pf.h.grids[0]._id_offset
        if self.comm.size > len(pf.h.grids):
            mylog.info('Parallel rendering requires that the number of \n \
            grids in the dataset is greater or equal to the number of \n \
            processors.  Reduce number of processors.')
            raise(KeyError)
        if fields is None: fields = ["Density"]
        self.no_ghost = no_ghost
        reduction_needed = {'domain':False,'depth':True,'breadth':True}
        self.tree_type = 'domain' # Hard code for now tree_type
        self.reduce_tree=reduction_needed[self.tree_type]
        self.bricks_loaded = False
        self.bricks = []
        self.brick_dimensions = []
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
            self.domain_left_edge = na.array(le)

        if re is None:
            self.domain_right_edge = pf.domain_right_edge
        else:
            self.domain_right_edge = na.array(re)

        self.domain_left_edge = na.clip(self.domain_left_edge,pf.domain_left_edge, pf.domain_right_edge)
        self.domain_right_edge = na.clip(self.domain_right_edge,pf.domain_left_edge, pf.domain_right_edge)

        levels = pf.hierarchy.get_levels()
        root_grids = levels.next()
        covering_grids = root_grids
        vol_needed = na.prod(self.domain_right_edge-self.domain_left_edge)

        for i in range(self.pf.hierarchy.max_level):
            root_l_data = na.clip(na.array([grid.LeftEdge for grid in root_grids]),self.domain_left_edge, self.domain_right_edge)
            root_r_data = na.clip(na.array([grid.RightEdge for grid in root_grids]),self.domain_left_edge, self.domain_right_edge)
            
            vol = na.prod(root_r_data-root_l_data,axis=1).sum()
            if vol >= vol_needed:
                covering_grids = root_grids
                root_grids = levels.next()
            else:
                break
            
        root_grids = covering_grids
        
        rgdds = root_grids[0].dds
        self.domain_left_edge = ((self.domain_left_edge)/rgdds).astype('int64')*rgdds
        self.domain_right_edge = (((self.domain_right_edge)/rgdds).astype('int64')+1)*rgdds

        self.domain_left_edge = na.clip(self.domain_left_edge,pf.domain_left_edge, pf.domain_right_edge)
        self.domain_right_edge = na.clip(self.domain_right_edge,pf.domain_left_edge, pf.domain_right_edge)
        
        self.my_l_corner = self.domain_left_edge
        self.my_r_corner = self.domain_right_edge

        #mylog.info('Making kd tree from le %s to %s'% (self.domain_left_edge, self.domain_right_edge))
        
        root_l_data = na.array([grid.LeftEdge for grid in root_grids])
        root_r_data = na.array([grid.RightEdge for grid in root_grids])
        root_we_want = na.all(root_l_data < self.my_r_corner,axis=1)*\
                       na.all(root_r_data > self.my_l_corner,axis=1)
        
        root_grids = root_grids[root_we_want]

        # Build the kd-Tree
        t1 = time()
        self._build(root_grids, None, self.domain_left_edge, self.domain_right_edge)
        t2 = time()
        mylog.debug('It took %e seconds to build AMRKDTree.tree' % (t2-t1))
        
        self._build_tree_dict()

        # If the full amr kD-tree is requested, merge the results from
        # the parallel build.
        if merge_trees and self.comm.size > 1:
            self.join_parallel_trees()            
            self.my_l_corner = self.domain_left_edge
            self.my_r_corner = self.domain_right_edge
        
        # Initialize the kd leafs:
        self.initialize_leafs()
        
        # Add properties to leafs/nodes
        self.total_cost = self.count_cost()

        # Calculate the total volume spanned by the tree
        self.volume = self.count_volume()
        #mylog.debug('Cost is %d' % self.total_cost)
        mylog.debug('Volume is %e' % self.volume) 

        self.current_saved_grids = []
        self.current_vcds = []


    def _build_tree_dict(self):
        self.tree_dict = {}
        for node in self.depth_traverse():
            self.tree_dict[node.id] = node

    def _overlap_check(self, le, re, brick, periodic=True):
        r"""Given a left and right edges along with a brick, tests overlap of any
        cells in the brick

        Parameters
        ----------
        le: array_like
            The left edge of the region being searched for overlap.
        re: array_like
            The right edge of the region being searched for overlap.
        periodic: boolean, optional
            Specifies whether search should include periodicity.  Default:True

        Returns
        ----------
        boolean: True if overlap is found, False otherwise.
        
        """
        if (le[0] < brick.r_corner[0]) and (re[0] > brick.l_corner[0]) and \
               (le[1] < brick.r_corner[1]) and (re[1] > brick.l_corner[1]) and \
               (le[2] < brick.r_corner[2]) and (re[2] > brick.l_corner[2]):
            return True

        if periodic:
            myle = deepcopy(le)
            myre = deepcopy(re)
            w = self.pf.domain_right_edge-self.pf.domain_left_edge
            for i in range(3):
                if myle[i] < self.pf.domain_left_edge[i]:
                    myle[i] += w[i]
                    myre[i] += w[i]
                if myre[i] > self.pf.domain_right_edge[i]:
                    myle[i] -= w[i]
                    myre[i] -= w[i]
                    
            if (myle[0] < brick.r_corner[0]) and (myre[0] > brick.l_corner[0]) and \
                   (myle[1] < brick.r_corner[1]) and (myre[1] > brick.l_corner[1]) and \
                   (myle[2] < brick.r_corner[2]) and (myre[2] > brick.l_corner[2]):
                return True
                
        return False

    def get_all_neighbor_bricks(self, node, le=None, re=None, periodic=True, add_to_brick_obj=False):
        r"""Given a brick_id, finds all other bricks that share a face, edge, or
        vertex.  Alternatively, will find all neighbors to an
        arbitrary rectangular specified by left and right edges.

        Parameters
        ----------
        brick_id: int
            ID of the brick in question.
        le: array_like, optional
            The left edge of an arbitrarily specified rectangular solid
        re: array_like, optional
            The right edge of an arbitrarily specified rectangular solid
        periodic: boolean, optional
            Specifies whether search should include periodicity.  Default:True
        iterator: boolean, optional
            If true, will yield the brick ids instead of return a list

        Returns
        ----------
        neighbors: list
           A list of all neighbor brick ids.
        
        """
        neighbors = []
        dx = self.pf.h.get_smallest_dx()
        if le is None:
            le = node.l_corner - dx
        if re is None:
            re = node.r_corner + dx

        nodes_to_check = [self.tree]
        while len(nodes_to_check) > 0:
            thisnode = nodes_to_check.pop(0)
            if thisnode.grid is None:
                if self._overlap_check(le,re,thisnode.left_child,periodic=periodic):
                    nodes_to_check.append(thisnode.left_child)
                if self._overlap_check(le,re,thisnode.right_child,periodic=periodic):
                    nodes_to_check.append(thisnode.right_child)
            else:
                neighbors.append(thisnode)

        if add_to_brick_obj:
            self.tree.neighbor_bricks=neighbors
        return neighbors

    def get_all_neighbor_grids(self, brick, le=None, re=None, periodic=True):
        r"""Given a brick_id, finds all other grids that share a face, edge, or
        vertex.  Alternatively, will find all neighbors to an
        arbitrary rectangular specified by left and right edges.

        Parameters
        ----------
        brick_id: int
            ID of the brick in question.
        le: array_like, optional
            The left edge of an arbitrarily specified rectangular solid
        re: array_like, optional
            The right edge of an arbitrarily specified rectangular solid
        periodic: boolean, optional
            Specifies whether search should include periodicity.  Default:True
        iterator: boolean, optional
            If true, will yield the grid ids instead of return a list

        Returns
        ----------
        neighbors: list
           A list of all neighbor grid ids.
        
        """
        grids = [brick.grid for brick in self.get_all_neighbor_bricks(
            brick, le=le, re=re, periodic=periodic)]
        return grids

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
        position = na.array(position)
        grid = self.locate_brick(position).grid
        ci = ((position-grid.LeftEdge)/grid.dds).astype('int64')
        return self.locate_neighbors(grid,ci)

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
        ci = na.array(ci)
        center_dds = grid.dds
        position = grid.LeftEdge + (na.array(ci)+0.5)*grid.dds
        grids = na.empty(26, dtype='object')
        cis = na.empty([26,3], dtype='int64')
        offs = 0.5*(center_dds + self.sdx)

        new_cis = ci + steps
        in_grid = na.all((new_cis >=0)*
                         (new_cis < grid.ActiveDimensions),axis=1)
        new_positions = position + steps*offs
        grids[in_grid] = grid
                
        get_them = na.argwhere(in_grid != True).ravel()
        cis[in_grid] = new_cis[in_grid]

        if (in_grid != True).sum()>0:
            grids[in_grid != True] = \
                [self.locate_brick(new_positions[i]).grid for i in get_them]
            cis[in_grid != True] = \
                [(new_positions[i]-grids[i].LeftEdge)/
                 grids[i].dds for i in get_them]
        cis = [tuple(ci) for ci in cis]
        return grids, cis

    def locate_brick(self, position):
        r"""Given a position, find the node that contains it.

        Will modify the position to account for periodicity.
        
        Parameters
        ----------
        pos: array_like
            Position being queried

        Returns
        ----------
        node: MasterNode()
            AMRKDTree node that contains position.
        
        """
        node = self.tree
        w = self.pf.domain_right_edge-self.pf.domain_left_edge
        for i in range(3):
            if position[i] < self.pf.domain_left_edge[i]:
                position[i] += w[i]
            if position[i] > self.pf.domain_right_edge[i]:
                position[i] -= w[i]
        while True:
            if node.grid is not None:
                return node
            else:
                if position[node.split_ax] < node.split_pos:
                    node = node.left_child
                else:
                    node = node.right_child
        return node
    
    def initialize_source(self):
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

        for current_node in self.depth_traverse():
            if current_node.grid is not None:
                if current_node.grid in current_saved_grids:
                    dds = current_vcds[current_saved_grids.index(current_node.grid)]
                else:
                    dds = []
                    for i,field in enumerate(self.fields):
                        vcd = current_node.grid.get_vertex_centered_data(field,smoothed=True,no_ghost=self.no_ghost).astype('float64')
                        if self.log_fields[i]: vcd = na.log10(vcd)
                        dds.append(vcd)
                    current_saved_grids.append(current_node.grid)
                    current_vcds.append(dds)

                data = [d[current_node.li[0]:current_node.ri[0]+1,
                          current_node.li[1]:current_node.ri[1]+1,
                          current_node.li[2]:current_node.ri[2]+1].copy() for d in dds]
                
                if na.any(current_node.r_corner-current_node.l_corner == 0):
                    current_node.brick = None
                else:
                    current_node.brick = PartitionedGrid(current_node.grid.id, data,
                                                         current_node.l_corner.copy(), 
                                                         current_node.r_corner.copy(), 
                                                         current_node.dims.astype('int64'))
                self.bricks.append(current_node.brick)
                self.brick_dimensions.append(current_node.dims)
        self.bricks = na.array(self.bricks)
        self.brick_dimensions = na.array(self.brick_dimensions)
        del current_saved_grids, current_vcds
        self.bricks_loaded = True

    def get_brick_data(self,current_node):
        if current_node.brick is not None:
            return 

        if current_node.grid in self.current_saved_grids:
            dds = self.current_vcds[self.current_saved_grids.index(current_node.grid)]
        else:
            dds = []
            for i,field in enumerate(self.fields):
                vcd = current_node.grid.get_vertex_centered_data(field,smoothed=True,no_ghost=self.no_ghost).astype('float64')
                if self.log_fields[i]: vcd = na.log10(vcd)
                dds.append(vcd)
                self.current_saved_grids.append(current_node.grid)
                self.current_vcds.append(dds)
                
        data = [d[current_node.li[0]:current_node.ri[0]+1,
                  current_node.li[1]:current_node.ri[1]+1,
                  current_node.li[2]:current_node.ri[2]+1].copy() for d in dds]

        current_node.brick = PartitionedGrid(current_node.grid.id, data,
                                             current_node.l_corner.copy(), 
                                             current_node.r_corner.copy(), 
                                             current_node.dims.astype('int64'))

        return
        
    def set_leaf_props(self,thisnode):
        r"""Given a leaf, gathers grid, indices, dimensions, and cost properties.

        Parameters
        ----------
        None

        Returns
        ----------
        None
        
        """
        thisnode.grid = self.pf.hierarchy.grids[thisnode.grid - self._id_offset]
        
        dds = thisnode.grid.dds
        gle = thisnode.grid.LeftEdge
        gre = thisnode.grid.RightEdge
        thisnode.li = na.rint((thisnode.l_corner-gle)/dds).astype('int32')
        thisnode.ri = na.rint((thisnode.r_corner-gle)/dds).astype('int32')
        thisnode.dims = (thisnode.ri - thisnode.li).astype('int32')
        # Here the cost is actually inversely proportional to 4**Level (empirical)
        #thisnode.cost = (na.prod(thisnode.dims)/4.**thisnode.grid.Level).astype('int64')
        thisnode.cost = 1.0
        # Here is the old way
        # thisnode.cost = na.prod(thisnode.dims).astype('int64')

    def initialize_leafs(self):
        for node in self.depth_traverse():
            if node.grid is not None:
                self.set_leaf_props(node)

    def join_parallel_trees(self):
        self.trim_references()
        self.merge_trees()
        self.rebuild_references()
                
    def trim_references(self):
        par_tree_depth = long(na.log2(self.comm.size))
        for i in range(2**self.comm.size):
            if ((i + 1)>>par_tree_depth) == 1:
                # There are self.comm.size nodes that meet this criteria
                if (i+1-self.comm.size) != self.comm.rank:
                    self.tree_dict.pop(i)
                    continue
        for node in self.tree_dict.itervalues():
            del node.parent, node.left_child, node.right_child
            try:
                del node.grids
            except:
                pass
            if not na.isreal(node.grid):
                node.grid = node.grid.id
        if self.tree_dict[0].split_pos is None:
            self.tree_dict.pop(0)
    def merge_trees(self):
        self.tree_dict = self.comm.par_combine_object(self.tree_dict,
                            datatype = "dict", op = "join")

    def rebuild_references(self):
        self.tree = self.tree_dict[0]
        self.tree.parent = None
        for node in self.depth_traverse():
            try:
                node.parent = self.tree_dict[_parent_id(node.id)]
            except:
                node.parent = None
            try:
                node.left_child = self.tree_dict[_lchild_id(node.id)]
            except:
                node.left_child = None
            try:
                node.right_child = self.tree_dict[_rchild_id(node.id)]
            except:
                node.right_child = None

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

        for node in self.depth_traverse():
            if node.grid is None:
                try:
                    node.cost = node.left_child.cost 
                except:
                    node.cost = 0
                try: node.cost += node.right_child.cost
                except:
                    pass
        return self.tree.cost

    def depth_traverse(self):
        '''
        Yields a depth-first traversal of the kd tree always going to
        the left child before the right.
        '''
        current = self.tree
        previous = None
        while current is not None:
            yield current
            current, previous = self.step_depth(current, previous)

    def step_depth(self, current, previous):
        '''
        Takes a single step in the depth-first traversal
        '''
        if current.grid is not None: # At a leaf, move back up
            previous = current
            # mylog.debug('moving up from leaf')
            current = current.parent
            
        elif current.parent is previous: # Moving down, go left first
            previous = current
            if current.left_child is not None:
                # mylog.debug('moving down to left child')
                current = current.left_child
            elif current.right_child is not None:
                # mylog.debug('no left, moving down to right child')
                current = current.right_child
            else:
                # mylog.debug('no left or right, moving to parent')
                current = current.parent
                
        elif current.left_child is previous: # Moving up from left, go right 
            previous = current
            if current.right_child is not None:
                # mylog.debug('moving down to right child')
                current = current.right_child
            else:
                # mylog.debug('no right, moving to parent')
                current = current.parent

        elif current.right_child is previous: # Moving up from right child, move up
            previous = current
            # mylog.debug('moving up from right child')
            current = current.parent
            
        return current, previous
    
    def viewpoint_traverse(self, viewpoint):
        '''
        Yields a viewpoint dependent traversal of the kd-tree.  Starts
        with nodes furthest away from viewpoint.
        '''
        
        current = self.tree
        previous = None
        # mylog.debug('Starting with %s %s'%(current, previous))
        while current is not None:
            yield current
            current, previous = self.step_viewpoint(current, previous, viewpoint)

    def step_viewpoint(self, current, previous, viewpoint):
        '''
        Takes a single step in the viewpoint based traversal.  Always
        goes to the node furthest away from viewpoint first.
        '''
        if current.grid is not None: # At a leaf, move back up
            previous = current
            current = current.parent
        elif current.split_ax is None: # This is a dead node
            previous = current
            current = current.parent

        elif current.parent is previous: # Moving down
            previous = current
            if viewpoint[current.split_ax] <= current.split_pos:
                if current.right_child is not None:
                    current = current.right_child
                else:
                    previous = current.right_child
            else:
                if current.left_child is not None:
                    current = current.left_child
                else:
                    previous = current.left_child
                
        elif current.right_child is previous: # Moving up from right 
            previous = current
            if viewpoint[current.split_ax] <= current.split_pos:
                if current.left_child is not None:
                    current = current.left_child
                else:
                    current = current.parent
            else:
                current = current.parent
                    
        elif current.left_child is previous: # Moving up from left child
            previous = current
            if viewpoint[current.split_ax] > current.split_pos:
                if current.right_child is not None:
                    current = current.right_child
                else:
                    current = current.parent
            else:
                current = current.parent
        
        return current, previous
                
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
        for node in self.depth_traverse():
            if node.grid is not None:
                v += na.prod(node.r_corner - node.l_corner)
        return v

    def count_cells(self):
        r"""Calculates the numbers of cells of the kd-Tree

        Parameters
        ----------
        None

        Returns
        ----------
        Total volume of the tree.
        
        """
        c = na.int64(0)
        for node in self.depth_traverse():
            if node.grid is not None:
                c += na.prod(node.ri - node.li).astype('int64')
        return c

    def _build(self, grids, parent, l_corner, r_corner):
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
        self.tree = MasterNode()

        head_node = self.tree
        previous_node = None
        current_node = self.tree
        current_node.grids = grids
        current_node.l_corner = l_corner
        current_node.r_corner = r_corner
        # current_node.owner = self.comm.rank
        current_node.id = 0
        par_tree_depth = int(na.log2(self.comm.size))
        anprocs = 2**par_tree_depth

        volume_partitioned = 0.0
        pbar = get_pbar("Building kd-Tree",
                na.prod(self.domain_right_edge-self.domain_left_edge))

        while current_node is not None:
            pbar.update(volume_partitioned)

            # If we don't have any grids, that means we are revisiting
            # a dividing node, and there is nothing to be done.
            try: ngrids = current_node.grids
            except:
                current_node, previous_node = self.step_depth(current_node, previous_node)
                continue

            # This is where all the domain decomposition occurs.  
            if ((current_node.id + 1)>>par_tree_depth) == 1:
                # There are anprocs nodes that meet this criteria
                if (current_node.id+1-anprocs) == self.comm.rank:
                    # I own this shared node
                    self.my_l_corner = current_node.l_corner
                    self.my_r_corner = current_node.r_corner
                else:
                    # This node belongs to someone else, move along
                    current_node, previous_node = self.step_depth(current_node, previous_node)
                    continue

            # If we are down to one grid, we are either in it or the parent grid
            if len(current_node.grids) == 1:
                thisgrid = current_node.grids[0]
                # If we are completely contained by that grid
                if (thisgrid.LeftEdge[0] <= current_node.l_corner[0]) and (thisgrid.RightEdge[0] >= current_node.r_corner[0]) and \
                   (thisgrid.LeftEdge[1] <= current_node.l_corner[1]) and (thisgrid.RightEdge[1] >= current_node.r_corner[1]) and \
                   (thisgrid.LeftEdge[2] <= current_node.l_corner[2]) and (thisgrid.RightEdge[2] >= current_node.r_corner[2]):
                    # Check if we have children and have not exceeded l_max
                    if len(thisgrid.Children) > 0 and thisgrid.Level < self.l_max:
                        # Get the children that are actually in the current volume
                        children = [child.id - self._id_offset for child in thisgrid.Children  
                                    if na.all(child.LeftEdge < current_node.r_corner) & 
                                    na.all(child.RightEdge > current_node.l_corner)]

                        # If we have children, get all the new grids, and keep building the tree
                        if len(children) > 0:
                            current_node.grids = self.pf.hierarchy.grids[na.array(children,copy=False)]
                            current_node.parent_grid = thisgrid
                            #print 'My single grid covers the rest of the volume, and I have children, about to iterate on them'
                            del children
                            continue

                    # Else make a leaf node (brick container)
                    #print 'My single grid covers the rest of the volume, and I have no children', thisgrid
                    set_leaf(current_node, thisgrid, current_node.l_corner, current_node.r_corner)
                    volume_partitioned += na.prod(current_node.r_corner-current_node.l_corner)
                    # print 'My single grid covers the rest of the volume, and I have no children'
                    current_node, previous_node = self.step_depth(current_node, previous_node)
                    continue

            # If we don't have any grids, this volume belongs to the parent        
            if len(current_node.grids) == 0:
                #print 'This volume does not have a child grid, so it belongs to my parent!'
                set_leaf(current_node, current_node.parent_grid, current_node.l_corner, current_node.r_corner)
                current_node, previous_node = self.step_depth(current_node, previous_node)
                continue

            # If we've made it this far, time to build a dividing node
            # print 'Building dividing node'
            # Continue if building failed
            if self._build_dividing_node(current_node): continue

            # Step to the nest node in a depth-first traversal.
            current_node, previous_node = self.step_depth(current_node, previous_node)
        
        pbar.finish()


    def _get_choices(self, current_node):
        '''
        Given a node, finds all the choices for the next dividing plane.  
        '''
        # For some reason doing dim 0 separately is slightly faster.
        # This could be rewritten to all be in the loop below.

        data = na.array([(child.LeftEdge, child.RightEdge) for child in current_node.grids],copy=False)
        best_dim, split, less_ids, greater_ids = \
            kdtree_get_choices(data, current_node.l_corner, current_node.r_corner)
        return data[:,:,best_dim], best_dim, split, less_ids, greater_ids

    def _build_dividing_node(self, current_node):
        '''
        Makes the current node a dividing node, and initializes the
        left and right children.
        '''

        data = na.array([(child.LeftEdge, child.RightEdge) for child in current_node.grids],copy=False)
        best_dim, split, less_ids, greater_ids = \
            kdtree_get_choices(data, current_node.l_corner, current_node.r_corner)

        del data

        # Here we break out if no unique grids were found. In this case, there
        # are likely overlapping grids, and we assume that the first grid takes
        # precedence.  This is fragile.
        if best_dim == -1:
            current_node.grids = [current_node.grids[0]]
            return 1

        current_node.split_ax = best_dim
        current_node.split_pos = split
        #less_ids0 = (data[:,0] < split)
        #greater_ids0 = (split < data[:,1])
        #assert(na.all(less_ids0 == less_ids))
        #assert(na.all(greater_ids0 == greater_ids))

        current_node.left_child = MasterNode(my_id=_lchild_id(current_node.id),
                                             parent=current_node,
                                             parent_grid=current_node.parent_grid,
                                             grids=current_node.grids[less_ids],
                                             l_corner=current_node.l_corner,
                                             r_corner=corner_bounds(best_dim, split, current_right=current_node.r_corner))

        current_node.right_child = MasterNode(my_id=_rchild_id(current_node.id),
                                              parent=current_node,
                                              parent_grid=current_node.parent_grid,
                                              grids=current_node.grids[greater_ids],
                                              l_corner=corner_bounds(best_dim, split, current_left=current_node.l_corner),
                                              r_corner=current_node.r_corner)

        # You have to at least delete the .grids object for the kd
        # build to work.  The other deletions are just to save memory.
        del current_node.grids, current_node.parent_grid, current_node.brick,\
            current_node.li, current_node.ri, current_node.dims

        return 0

    def traverse(self, back_center, front_center, image):
        r"""Traverses the kd-Tree, casting the partitioned grids from back to
            front.

        Given the back and front centers, and the image, ray-cast
        using the kd-Tree structure.

        Parameters
        ----------
        back_center: array_like
            Position of the back center from which to start moving forward.
        front_center: array_like
            Position of the front center to which the traversal progresses.
        image: na.array
            Image plane to contain resulting ray cast.

        Returns
        ----------
        None, but modifies the image array.

        See Also
        ----------
        yt.visualization.volume_rendering.camera

        """
        if self.tree is None: 
            print 'No KD Tree Exists'
            return
        self.image = image

        viewpoint = back_center 
        print 'Moving from front_center to back_center:',front_center, back_center

        for node in self.viewpoint_traverse(viewpoint):
            if node.grid is not None:
                if node.brick is not None:
                    yield node.brick

        mylog.debug('About to enter reduce, my image has a max of %e' %
                self.image.max())
        self.reduce_tree_images(self.tree, viewpoint)
        self.comm.barrier()

    def reduce_tree_images(self, tree, viewpoint, image=None):
        if image is not None:
            self.image = image
        rounds = int(na.log2(self.comm.size))
        anprocs = 2**rounds
        my_node = tree
        my_node_id = 0
        my_node.owner = 0
        path = na.binary_repr(anprocs+self.comm.rank)
        for i in range(rounds):
            try:
                my_node.left_child.owner = my_node.owner
                my_node.right_child.owner = my_node.owner + 2**(rounds-(i+1))
                if path[i+1] == '0':
                    my_node = my_node.left_child
                    my_node_id = my_node.id
                else:
                    my_node = my_node.right_child
                    my_node_id = my_node.id
            except:
                rounds = i-1
        for thisround in range(rounds,0,-1):
            #print self.comm.rank, 'my node', my_node_id
            parent = my_node.parent
            #print parent['split_ax'], parent['split_pos']
            if viewpoint[parent.split_ax] > parent.split_pos:
                front = parent.right_child
                back = parent.left_child
            else:
                front = parent.left_child
                back = parent.right_child
            print 'Combining', viewpoint, parent.split_ax, parent.split_pos
            print front.l_corner, front.r_corner
            print back.l_corner, back.r_corner

            # mylog.debug('front owner %i back owner %i parent owner %i'%( front.owner, back.owner, parent.owner))
            # Send the images around
            if front.owner == self.comm.rank:
                if front.owner == parent.owner:
                    mylog.debug( '%04i receiving image from %04i'%(self.comm.rank,back.owner))
                    arr2 = self.comm.recv_array(back.owner, tag=back.owner).reshape(
                        (self.image.shape[0],self.image.shape[1],self.image.shape[2]))
                    ta = 1.0 - na.sum(self.image,axis=2)
                    ta[ta<0.0] = 0.0
                    for i in range(3):
                        # This is the new way: alpha corresponds to opacity of a given
                        # slice.  Previously it was ill-defined, but represented some
                        # measure of emissivity.
                        self.image[:,:,i  ] = self.image[:,:,i] + ta*arr2[:,:,i]
                else:
                    mylog.debug('Reducing image.  You have %i rounds to go in this binary tree' % thisround)
                    mylog.debug('%04i sending my image to %04i with max %e'%(self.comm.rank,back.owner, self.image.max()))
                    self.comm.send_array(self.image.ravel(), back.owner, tag=self.comm.rank)

            if back.owner == self.comm.rank:
                if front.owner == parent.owner:
                    mylog.debug('%04i sending my image to %04i with max %e'%(self.comm.rank, front.owner,
                                self.image.max()))
                    self.comm.send_array(self.image.ravel(), front.owner, tag=self.comm.rank)
                else:
                    mylog.debug('Reducing image.  You have %i rounds to go in this binary tree' % thisround)
                    mylog.debug('%04i receiving image from %04i'%(self.comm.rank,front.owner))
                    arr2 = self.comm.recv_array(front.owner, tag=front.owner).reshape(
                        (self.image.shape[0],self.image.shape[1],self.image.shape[2]))
                    #ta = na.exp(-na.sum(arr2,axis=2))
                    ta = 1.0 - na.sum(arr2, axis=2)
                    ta[ta<0.0] = 0.0
                    for i in range(3):
                        # This is the new way: alpha corresponds to opacity of a given
                        # slice.  Previously it was ill-defined, but represented some
                        # measure of emissivity.
                        self.image[:,:,i  ] = ta*self.image[:,:,i  ] + arr2[:,:,i]

            # Set parent owner to back owner
            # my_node = (my_node-1)>>1
            if self.comm.rank == my_node.parent.owner: 
                my_node = my_node.parent
            else:
                break

    def store_kd_bricks(self, fn=None):
        if fn is None:
            fn = '%s_kd_bricks.h5'%self.pf
        if self.comm.rank != 0:
            self.comm.recv_array(self.comm.rank-1, tag=self.comm.rank-1)
        f = h5py.File(fn,"a")
        for node in self.depth_traverse():
            i = node.id
            if node.grid is not None:
                for fi,field in enumerate(self.fields):
                    try:
                        f.create_dataset("/brick_%s_%s" % (hex(i),field),
                                         data = node.brick.my_data[fi].astype('float64'))
                    except:
                        pass
        f.close()
        if self.comm.rank != (self.comm.size-1):
            self.comm.send_array([0],self.comm.rank+1, tag=self.comm.rank)
        
    def load_kd_bricks(self,fn=None):
        if fn is None:
            fn = '%s_kd_bricks.h5' % self.pf
        if self.comm.rank != 0:
            self.comm.recv_array(self.comm.rank-1, tag=self.comm.rank-1)
        try:
            f = h5py.File(fn,"r")
            for node in self.depth_traverse():
                i = node.id
                if node.grid is not None:
                    data = [f["brick_%s_%s" %
                              (hex(i), field)][:].astype('float64') for field in self.fields]
                    node.brick = PartitionedGrid(node.grid.id, data,
                                                 node.l_corner.copy(), 
                                                 node.r_corner.copy(), 
                                                 node.dims.astype('int64'))
                    
                    self.bricks.append(node.brick)
                    self.brick_dimensions.append(node.dims)

            self.bricks = na.array(self.bricks)
            self.brick_dimensions = na.array(self.brick_dimensions)

            self.bricks_loaded=True
            f.close()
        except:
            pass
        if self.comm.rank != (self.comm.size-1):
            self.comm.send_array([0],self.comm.rank+1, tag=self.comm.rank)

    def load_tree(self,fn):
        raise NotImplementedError()
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
        raise NotImplementedError()
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
        
    def corners_to_line(self,lc, rc):
        x = na.array([ lc[0], lc[0], lc[0], lc[0], lc[0],
                       rc[0], rc[0], rc[0], rc[0], rc[0],
                       rc[0], lc[0], lc[0], rc[0],
                       rc[0], lc[0], lc[0] ])
        
        y = na.array([ lc[1], lc[1], rc[1], rc[1], lc[1],
                       lc[1], lc[1], rc[1], rc[1], lc[1],
                       lc[1], lc[1], rc[1], rc[1],
                       rc[1], rc[1], lc[1] ])
        
        z = na.array([ lc[2], rc[2], rc[2], lc[2], lc[2],
                       lc[2], rc[2], rc[2], lc[2], lc[2],
                       rc[2], rc[2], rc[2], rc[2],
                       lc[2], lc[2], lc[2] ])
        return x,y,z
