"""
Tests for GridTree

Author: John ZuHone <jzuhone@gmail.com>
Affiliation: NASA/Goddard Space Flight Center
Homepage: http://yt-project.org/
License:
Copyright (C) 2012 John ZuHone.  All Rights Reserved.

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
        grid["Density"] = \
            np.random.random(grid["dimensions"]) * 2 ** grid["level"]
    test_pf = load_amr_grids(grid_data, [16, 16, 16], 1.0)


def test_grid_tree():
    """Main test suite for GridTree"""
    grid_tree = test_pf.h.get_grid_tree()
    indices, levels, nchild, children = grid_tree.return_tree_info()

    grid_levels = [grid.Level for grid in test_pf.h.grids]

    grid_indices = [grid.id - grid._id_offset for grid in test_pf.h.grids]
    grid_nchild = [len(grid.Children) for grid in test_pf.h.grids]

    yield assert_equal, levels, grid_levels
    yield assert_equal, indices, grid_indices
    yield assert_equal, nchild, grid_nchild

    for i, grid in enumerate(test_pf.h.grids):
        if grid_nchild[i] > 0:
            grid_children = np.array([child.id - child._id_offset
                                      for child in grid.Children])
            yield assert_equal, grid_children, children[i]


def test_find_points():
    """Main test suite for MatchPoints"""
    num_points = 100
    randx = np.random.uniform(low=test_pf.domain_left_edge[0],
                              high=test_pf.domain_right_edge[0],
                              size=num_points)
    randy = np.random.uniform(low=test_pf.domain_left_edge[1],
                              high=test_pf.domain_right_edge[1],
                              size=num_points)
    randz = np.random.uniform(low=test_pf.domain_left_edge[2],
                              high=test_pf.domain_right_edge[2],
                              size=num_points)

    point_grids, point_grid_inds = test_pf.h.find_points(randx, randy, randz)

    grid_inds = np.zeros((num_points), dtype='int64')

    for ind, ixx, iyy, izz in zip(range(num_points), randx, randy, randz):

        pt_level = -1

        for grid in test_pf.h.grids:

            if grid.is_in_grid(ixx, iyy, izz):

                if grid.Level > pt_level:
                    pt_level = grid.Level
                    grid_inds[ind] = grid.id - grid._id_offset

    yield assert_equal, point_grid_inds, grid_inds

    # Test wheter find_points works for lists
    point_grids, point_grid_inds = test_pf.h.find_points(randx.tolist(),
                                                         randy.tolist(),
                                                         randz.tolist())
    yield assert_equal, point_grid_inds, grid_inds

    # Test if find_points works for scalar
    ind = random.randint(0, num_points - 1)
    point_grids, point_grid_inds = test_pf.h.find_points(randx[ind],
                                                         randy[ind],
                                                         randz[ind])
    yield assert_equal, point_grid_inds, grid_inds[ind]

    # Test if find_points fails properly for non equal indices' array sizes
    yield assert_raises, AssertionError, test_pf.h.find_points, \
        [0], 1.0, [2, 3]
