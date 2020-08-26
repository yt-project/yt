import tempfile
import time

import numpy as np
from nose.tools import assert_raises

import yt.utilities.lib.cykdtree as cykdtree
from yt.utilities.lib.cykdtree.tests import (
    make_points,
    make_points_neighbors,
    parametrize,
)


@parametrize(
    npts=100, ndim=(2, 3), periodic=(False, True), use_sliding_midpoint=(False, True)
)
def test_PyKDTree(npts=100, ndim=2, periodic=False, use_sliding_midpoint=False):
    pts, le, re, ls = make_points(npts, ndim)
    cykdtree.PyKDTree(
        pts,
        le,
        re,
        leafsize=ls,
        periodic=periodic,
        use_sliding_midpoint=use_sliding_midpoint,
    )


def test_PyKDTree_errors():
    pts, le, re, ls = make_points(100, 2)
    assert_raises(ValueError, cykdtree.PyKDTree, pts, le, re, leafsize=1)


@parametrize(npts=100, ndim=(2, 3), periodic=(False, True))
def test_search(npts=100, ndim=2, periodic=False):
    pts, le, re, ls = make_points(npts, ndim)
    tree = cykdtree.PyKDTree(pts, le, re, leafsize=ls, periodic=periodic)
    pos_list = [le, (le + re) / 2.0]
    if periodic:
        pos_list.append(re)
    for pos in pos_list:
        leaf = tree.get(pos)
        leaf.neighbors


@parametrize(npts=100, ndim=(2, 3))
def test_search_errors(npts=100, ndim=2):
    pts, le, re, ls = make_points(npts, ndim)
    tree = cykdtree.PyKDTree(pts, le, re, leafsize=ls)
    assert_raises(ValueError, tree.get, re)


@parametrize(periodic=(False, True))
def test_neighbors(periodic=False):
    pts, le, re, ls, left_neighbors, right_neighbors = make_points_neighbors(
        periodic=periodic
    )
    tree = cykdtree.PyKDTree(pts, le, re, leafsize=ls, periodic=periodic)
    for leaf in tree.leaves:
        out_str = str(leaf.id)
        try:
            for d in range(tree.ndim):
                out_str += "\nleft:  {} {} {}".format(
                    d, leaf.left_neighbors[d], left_neighbors[d][leaf.id]
                )
                assert len(left_neighbors[d][leaf.id]) == len(leaf.left_neighbors[d])
                for i in range(len(leaf.left_neighbors[d])):
                    assert left_neighbors[d][leaf.id][i] == leaf.left_neighbors[d][i]
                out_str += "\nright: {} {} {}".format(
                    d, leaf.right_neighbors[d], right_neighbors[d][leaf.id]
                )
                assert len(right_neighbors[d][leaf.id]) == len(leaf.right_neighbors[d])
                for i in range(len(leaf.right_neighbors[d])):
                    assert right_neighbors[d][leaf.id][i] == leaf.right_neighbors[d][i]
        except Exception as e:
            for leaf in tree.leaves:
                print(leaf.id, leaf.left_edge, leaf.right_edge)
            print(out_str)
            raise e


@parametrize(npts=100, ndim=(2, 3), periodic=(False, True))
def test_get_neighbor_ids(npts=100, ndim=2, periodic=False):
    pts, le, re, ls = make_points(npts, ndim)
    tree = cykdtree.PyKDTree(pts, le, re, leafsize=ls, periodic=periodic)
    pos_list = [le, (le + re) / 2.0]
    if periodic:
        pos_list.append(re)
    for pos in pos_list:
        tree.get_neighbor_ids(pos)


def time_tree_construction(Ntime, LStime, ndim=2):
    pts, le, re, ls = make_points(Ntime, ndim, leafsize=LStime)
    t0 = time.time()
    cykdtree.PyKDTree(pts, le, re, leafsize=LStime)
    t1 = time.time()
    print(f"{Ntime} {ndim}D points, leafsize {LStime}: took {t1 - t0} s")


def time_neighbor_search(Ntime, LStime, ndim=2):
    pts, le, re, ls = make_points(Ntime, ndim, leafsize=LStime)
    tree = cykdtree.PyKDTree(pts, le, re, leafsize=LStime)
    t0 = time.time()
    tree.get_neighbor_ids(0.5 * np.ones(tree.ndim, "double"))
    t1 = time.time()
    print(f"{Ntime} {ndim}D points, leafsize {LStime}: took {t1 - t0} s")


def test_save_load():
    for periodic in (True, False):
        for ndim in range(1, 5):
            pts, le, re, ls = make_points(100, ndim)
            tree = cykdtree.PyKDTree(
                pts, le, re, leafsize=ls, periodic=periodic, data_version=ndim + 12
            )
            with tempfile.NamedTemporaryFile(delete=False) as tf:
                tree.save(tf.name)
                restore_tree = cykdtree.PyKDTree.from_file(tf.name)
                tree.assert_equal(restore_tree)
