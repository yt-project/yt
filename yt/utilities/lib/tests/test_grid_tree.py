import numpy as np

from yt.testing import *
from yt.frontends.stream.api import load_amr_grids

def setup():

    global pf
    
    grid_data = [
        dict(left_edge = [0.0, 0.0, 0.0], right_edge = [1.0, 1.0, 1.],
             level = 0, dimensions = [16, 16, 16]),
        dict(left_edge = [0.25, 0.25, 0.25], right_edge = [0.75, 0.75, 0.75],
             level = 1, dimensions = [16, 16, 16]),
        dict(left_edge = [0.25, 0.25, 0.375], right_edge = [0.5, 0.5, 0.625],
             level = 2, dimensions = [16, 16, 16]),
        dict(left_edge = [0.5, 0.5, 0.375], right_edge = [0.75, 0.75, 0.625],
             level = 2, dimensions = [16, 16, 16]),
        dict(left_edge = [0.3125, 0.3125, 0.4375], right_edge = [0.4375, 0.4375, 0.5625],
             level = 3, dimensions = [16, 16, 16]),
        dict(left_edge = [0.5625, 0.5625, 0.4375], right_edge = [0.6875, 0.6875, 0.5625],
             level = 3, dimensions = [16, 16, 16])
        ]

    for g in grid_data: g["Density"] = np.random.random(g["dimensions"]) * 2**g["level"]
    pf = load_amr_grids(grid_data, [16, 16, 16], 1.0)

def test_grid_tree() :

    grid_tree = pf.h.get_grid_tree()
    indices, levels, nchild, children = grid_tree.return_tree_info()

    grid_levels = [grid.Level for grid in pf.h.grids]
    grid_indices = [grid.id-grid._id_offset for grid in pf.h.grids]
    grid_nchild = [len(grid.Children) for grid in pf.h.grids]

    print levels, grid_levels
    assert_equal(levels, grid_levels)
    assert_equal(indices, grid_indices)
    assert_equal(nchild, grid_nchild)

    for i, grid in enumerate(pf.h.grids) :
        if grid_nchild[i] > 0:
            grid_children = np.array([child.id-child._id_offset
                                      for child in grid.Children])
            assert_equal(grid_children, children[i])

def test_find_points() :
    
    num_points = 100

    x = np.random.uniform(low=pf.domain_left_edge[0],
                          high=pf.domain_right_edge[0], size=num_points)
    y = np.random.uniform(low=pf.domain_left_edge[1],
                          high=pf.domain_right_edge[1], size=num_points)
    z = np.random.uniform(low=pf.domain_left_edge[2],
                          high=pf.domain_right_edge[2], size=num_points)

    point_grids, point_grid_inds = pf.h.find_points(x,y,z)

    grid_inds = np.zeros((num_points), dtype='int64')

    for i, xx, yy, zz in zip(range(num_points), x, y, z) :

        pt_level = -1
        
        for grid in pf.h.grids:

            if grid.is_in_grid(xx, yy, zz) :
            
                if grid.Level > pt_level :
                    pt_level = grid.Level
                    grid_inds[i] = grid.id-grid._id_offset
                    
    assert_equal(point_grid_inds, grid_inds)
