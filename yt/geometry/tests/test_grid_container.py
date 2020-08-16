import random

import numpy as np

from yt.loaders import load_amr_grids
from yt.testing import assert_equal, assert_raises


def setup_test_ds():
    """Prepare setup specific environment"""
    grid_data = [
        dict(
            left_edge=[0.0, 0.0, 0.0],
            right_edge=[1.0, 1.0, 1.0],
            level=0,
            dimensions=[16, 16, 16],
        ),
        dict(
            left_edge=[0.25, 0.25, 0.25],
            right_edge=[0.75, 0.75, 0.75],
            level=1,
            dimensions=[16, 16, 16],
        ),
        dict(
            left_edge=[0.25, 0.25, 0.375],
            right_edge=[0.5, 0.5, 0.625],
            level=2,
            dimensions=[16, 16, 16],
        ),
        dict(
            left_edge=[0.5, 0.5, 0.375],
            right_edge=[0.75, 0.75, 0.625],
            level=2,
            dimensions=[16, 16, 16],
        ),
        dict(
            left_edge=[0.3125, 0.3125, 0.4375],
            right_edge=[0.4375, 0.4375, 0.5625],
            level=3,
            dimensions=[16, 16, 16],
        ),
        dict(
            left_edge=[0.5625, 0.5625, 0.4375],
            right_edge=[0.6875, 0.6875, 0.5625],
            level=3,
            dimensions=[16, 16, 16],
        ),
    ]

    for grid in grid_data:
        grid["density"] = (
            np.random.random(grid["dimensions"]) * 2 ** grid["level"],
            "g/cm**3",
        )
    return load_amr_grids(grid_data, [16, 16, 16])


def test_grid_tree():
    """Main test suite for GridTree"""
    test_ds = setup_test_ds()
    grid_tree = test_ds.index._get_grid_tree()
    indices, levels, nchild, children = grid_tree.return_tree_info()

    grid_levels = [grid.Level for grid in test_ds.index.grids]

    grid_indices = [grid.id - grid._id_offset for grid in test_ds.index.grids]
    grid_nchild = [len(grid.Children) for grid in test_ds.index.grids]

    assert_equal(levels, grid_levels)
    assert_equal(indices, grid_indices)
    assert_equal(nchild, grid_nchild)

    for i, grid in enumerate(test_ds.index.grids):
        if grid_nchild[i] > 0:
            grid_children = np.array(
                [child.id - child._id_offset for child in grid.Children]
            )
            assert_equal(grid_children, children[i])


def test_find_points():
    """Main test suite for MatchPoints"""
    num_points = 100
    test_ds = setup_test_ds()
    randx = np.random.uniform(
        low=test_ds.domain_left_edge[0],
        high=test_ds.domain_right_edge[0],
        size=num_points,
    )
    randy = np.random.uniform(
        low=test_ds.domain_left_edge[1],
        high=test_ds.domain_right_edge[1],
        size=num_points,
    )
    randz = np.random.uniform(
        low=test_ds.domain_left_edge[2],
        high=test_ds.domain_right_edge[2],
        size=num_points,
    )

    point_grids, point_grid_inds = test_ds.index._find_points(randx, randy, randz)

    grid_inds = np.zeros((num_points), dtype="int64")

    for ind, ixx, iyy, izz in zip(range(num_points), randx, randy, randz):

        pos = np.array([ixx, iyy, izz])
        pt_level = -1

        for grid in test_ds.index.grids:

            if (
                np.all(pos >= grid.LeftEdge)
                and np.all(pos <= grid.RightEdge)
                and grid.Level > pt_level
            ):
                pt_level = grid.Level
                grid_inds[ind] = grid.id - grid._id_offset

    assert_equal(point_grid_inds, grid_inds)

    # Test whether find_points works for lists
    point_grids, point_grid_inds = test_ds.index._find_points(
        randx.tolist(), randy.tolist(), randz.tolist()
    )
    assert_equal(point_grid_inds, grid_inds)

    # Test if find_points works for scalar
    ind = random.randint(0, num_points - 1)
    point_grids, point_grid_inds = test_ds.index._find_points(
        randx[ind], randy[ind], randz[ind]
    )
    assert_equal(point_grid_inds, grid_inds[ind])

    # Test if find_points fails properly for non equal indices' array sizes
    assert_raises(ValueError, test_ds.index._find_points, [0], 1.0, [2, 3])


def test_grid_arrays_view():
    ds = setup_test_ds()
    tree = ds.index._get_grid_tree()
    grid_arr = tree.grid_arrays
    assert_equal(grid_arr["left_edge"], ds.index.grid_left_edge)
    assert_equal(grid_arr["right_edge"], ds.index.grid_right_edge)
    assert_equal(grid_arr["dims"], ds.index.grid_dimensions)
    assert_equal(grid_arr["level"], ds.index.grid_levels[:, 0])
