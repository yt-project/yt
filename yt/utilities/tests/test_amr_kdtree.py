"""
Unit test the ARMKDTree in yt.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.amr_kdtree.api import AMRKDTree
from yt.utilities.lib.amr_kdtools import depth_traverse, \
        get_left_edge, get_right_edge
import yt.utilities.initial_conditions as ic
import yt.utilities.flagging_methods as fm
from yt.frontends.stream.api import load_uniform_grid, refine_amr
from yt.testing import assert_equal
import numpy as np


def test_amr_kdtree_coverage():
    return #TESTDISABLED
    domain_dims = (32, 32, 32)
    data = np.zeros(domain_dims) + 0.25
    fo = [ic.CoredSphere(0.05, 0.3, [0.7, 0.4, 0.75],
                         {"density": (0.25, 100.0)})]
    rc = [fm.flagging_method_registry["overdensity"](8.0)]
    ug = load_uniform_grid({"density": data}, domain_dims, 1.0)
    ds = refine_amr(ug, rc, fo, 5)

    kd = AMRKDTree(ds)

    volume = kd.count_volume()
    yield assert_equal, volume, \
        np.prod(ds.domain_right_edge - ds.domain_left_edge)

    cells = kd.count_cells()
    true_cells = ds.all_data().quantities['TotalQuantity']('Ones')[0]
    yield assert_equal, cells, true_cells

    # This largely reproduces the AMRKDTree.tree.check_tree() functionality
    tree_ok = True
    for node in depth_traverse(kd.tree.trunk):
        if node.grid is None:
            continue
        grid = ds.index.grids[node.grid - kd._id_offset]
        dds = grid.dds
        gle = grid.LeftEdge
        nle = get_left_edge(node)
        nre = get_right_edge(node)
        li = np.rint((nle-gle)/dds).astype('int32')
        ri = np.rint((nre-gle)/dds).astype('int32')
        dims = (ri - li).astype('int32')
        tree_ok *= np.all(grid.LeftEdge <= nle)
        tree_ok *= np.all(grid.RightEdge >= nre)
        tree_ok *= np.all(dims > 0)

    yield assert_equal, True, tree_ok
