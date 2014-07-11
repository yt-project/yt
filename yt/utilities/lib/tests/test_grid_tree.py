"""
Tests for GridTree



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import numpy as np
import random

from yt.testing import \
    assert_equal, assert_raises
from yt.frontends.stream.api import \
    load_amr_grids


def setup():
    """Prepare setup specific environment"""
    global test_pf

    grid_data = [
        dict(left_edge=[0.0, 0.0, 0.0], right_edge=[1.0, 1.0, 1.],
             level=0, dimensions=[16, 16, 16]),
        dict(left_edge=[0.25, 0.25, 0.25], right_edge=[0.75, 0.75, 0.75],
             level=1, dimensions=[16, 16, 16]),
        dict(left_edge=[0.25, 0.25, 0.375], right_edge=[0.5, 0.5, 0.625],
             level=2, dimensions=[16, 16, 16]),
        dict(left_edge=[0.5, 0.5, 0.375], right_edge=[0.75, 0.75, 0.625],
             level=2, dimensions=[16, 16, 16]),
        dict(left_edge=[0.3125, 0.3125, 0.4375],
             right_edge=[0.4375, 0.4375, 0.5625],
             level=3, dimensions=[16, 16, 16]),
        dict(left_edge=[0.5625, 0.5625, 0.4375],
             right_edge=[0.6875, 0.6875, 0.5625],
             level=3, dimensions=[16, 16, 16])
    ]

    for grid in grid_data:
        grid["density"] = \
            np.random.random(grid["dimensions"]) * 2 ** grid["level"]
    test_pf = load_amr_grids(grid_data, [16, 16, 16], 1.0)


def test_grid_tree():
    """Main test suite for GridTree"""
    grid_tree = test_pf.index.get_grid_tree()
    indices, levels, nchild, children = grid_tree.return_tree_info()

    grid_levels = [grid.Level for grid in test_pf.index.grids]

    grid_indices = [grid.id - grid._id_offset for grid in test_pf.index.grids]
    grid_nchild = [len(grid.Children) for grid in test_pf.index.grids]

    yield assert_equal, levels, grid_levels
    yield assert_equal, indices, grid_indices
    yield assert_equal, nchild, grid_nchild

    for i, grid in enumerate(test_pf.index.grids):
        if grid_nchild[i] > 0:
            grid_children = np.array([child.id - child._id_offset
                                      for child in grid.Children])
            yield assert_equal, grid_children, children[i]
