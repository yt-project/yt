import numpy as np
import pytest
from numpy.testing import assert_array_equal

from yt.testing import fake_amr_ds, fake_random_ds
from yt.utilities.amr_kdtree.api import AMRKDTree
from yt.utilities.lib.amr_kdtools import Node, viewpoint_node_ids

VIEWPOINTS = [
    (0.5, 0.5, 0.5),
    (0.1, 0.9, 0.3),
    (-10.0, 0.5, 0.5),
    (10.0, 0.5, 0.5),
    (0.5, -10.0, 0.5),
    (0.5, 10.0, 0.5),
    (0.5, 0.5, -10.0),
    (0.5, 0.5, 10.0),
    (-5.0, -5.0, -5.0),
    (5.0, 5.0, 5.0),
    (-5.0, 5.0, -5.0),
]


@pytest.fixture(scope="module", params=["amr", "uniform"])
def trunk(request):
    if request.param == "amr":
        ds = fake_amr_ds()
    else:
        ds = fake_random_ds(32, nprocs=8)
    trunk = AMRKDTree(ds).tree.trunk
    # node_ind is only populated by the contour finder, so give each leaf a
    # distinct value that differs from its node_id
    for i, node in enumerate(trunk.kd_traverse()):
        node.node_ind = 1000 + i
    return trunk


def _walk(trunk, viewpoint, size=None):
    n_leaves = len(list(trunk.kd_traverse()))
    if size is None:
        size = n_leaves
    node_ids = np.full(size, -1, dtype="int64")
    node_inds = np.full(size, -1, dtype="int64")
    n = viewpoint_node_ids(trunk, viewpoint, node_ids, node_inds)
    return n, node_ids, node_inds


@pytest.mark.parametrize("viewpoint", VIEWPOINTS)
def test_viewpoint_node_ids_matches_kd_traverse(trunk, viewpoint):
    expected = list(trunk.kd_traverse(viewpoint=np.array(viewpoint)))
    n, node_ids, node_inds = _walk(trunk, viewpoint)
    assert n == len(expected)
    assert_array_equal(node_ids, [node.node_id for node in expected])
    assert_array_equal(node_inds, [node.node_ind for node in expected])


@pytest.mark.parametrize("viewpoint", VIEWPOINTS)
def test_viewpoint_node_ids_visits_each_leaf_once(trunk, viewpoint):
    n, node_ids, _ = _walk(trunk, viewpoint)
    all_leaves = sorted(node.node_id for node in trunk.kd_traverse())
    assert n == len(all_leaves)
    assert_array_equal(np.sort(node_ids), all_leaves)


def test_viewpoint_node_ids_on_split_plane(trunk):
    # A viewpoint exactly on a split plane exercises the ``<=`` branch
    viewpoint = np.array([0.5, 0.5, 0.5])
    viewpoint[trunk.get_split_dim()] = trunk.get_split_pos()
    expected = [node.node_id for node in trunk.kd_traverse(viewpoint=viewpoint)]
    n, node_ids, _ = _walk(trunk, viewpoint)
    assert n == len(expected)
    assert_array_equal(node_ids, expected)


def test_viewpoint_node_ids_leaves_extra_buffer_untouched(trunk):
    n_leaves = len(list(trunk.kd_traverse()))
    n, node_ids, node_inds = _walk(trunk, (0.2, 0.4, 0.6), size=n_leaves + 5)
    assert n == n_leaves
    assert np.all(node_ids[:n_leaves] >= 0)
    assert np.all(node_inds[:n_leaves] >= 0)
    assert_array_equal(node_ids[n_leaves:], -1)
    assert_array_equal(node_inds[n_leaves:], -1)


@pytest.mark.parametrize(
    "viewpoint",
    [[0.1, 0.2, 0.3], (0.1, 0.2, 0.3), np.array([0.1, 0.2, 0.3])],
    ids=["list", "tuple", "ndarray"],
)
def test_viewpoint_node_ids_viewpoint_types(trunk, viewpoint):
    expected = [node.node_id for node in trunk.kd_traverse(viewpoint=viewpoint)]
    n, node_ids, _ = _walk(trunk, viewpoint)
    assert n == len(expected)
    assert_array_equal(node_ids, expected)


def _two_leaf_tree():
    le = np.zeros(3)
    re = np.ones(3)
    trunk = Node(None, None, None, le, re, -1, 1)
    trunk.create_split(0, 0.5)
    left_re = np.array([0.5, 1.0, 1.0])
    right_le = np.array([0.5, 0.0, 0.0])
    trunk.left = Node(trunk, None, None, le, left_re, 7, 2)
    trunk.right = Node(trunk, None, None, right_le, re, 8, 3)
    trunk.left.node_ind = 20
    trunk.right.node_ind = 30
    return trunk


@pytest.mark.parametrize(
    "x, expected_ids, expected_inds",
    [
        # Viewpoint on the left: the far (right) leaf comes first
        (0.25, [3, 2], [30, 20]),
        # Viewpoint on the right: the far (left) leaf comes first
        (0.75, [2, 3], [20, 30]),
        # On the split plane, ties go the same way as the left side
        (0.5, [3, 2], [30, 20]),
    ],
)
def test_viewpoint_node_ids_back_to_front(x, expected_ids, expected_inds):
    trunk = _two_leaf_tree()
    node_ids = np.empty(2, dtype="int64")
    node_inds = np.empty(2, dtype="int64")
    assert viewpoint_node_ids(trunk, (x, 0.5, 0.5), node_ids, node_inds) == 2
    assert_array_equal(node_ids, expected_ids)
    assert_array_equal(node_inds, expected_inds)


def test_viewpoint_node_ids_skips_gridless_leaves():
    trunk = _two_leaf_tree()
    trunk.left.grid = -1
    node_ids = np.full(2, -1, dtype="int64")
    node_inds = np.full(2, -1, dtype="int64")
    assert viewpoint_node_ids(trunk, (0.25, 0.5, 0.5), node_ids, node_inds) == 1
    assert_array_equal(node_ids, [3, -1])
    assert_array_equal(node_inds, [30, -1])


def test_viewpoint_node_ids_single_leaf():
    trunk = Node(None, None, None, np.zeros(3), np.ones(3), 4, 1)
    trunk.node_ind = 9
    node_ids = np.empty(1, dtype="int64")
    node_inds = np.empty(1, dtype="int64")
    assert viewpoint_node_ids(trunk, (0.5, 0.5, 0.5), node_ids, node_inds) == 1
    assert_array_equal(node_ids, [1])
    assert_array_equal(node_inds, [9])


def test_viewpoint_node_ids_rejects_wrong_dtype():
    trunk = _two_leaf_tree()
    with pytest.raises(ValueError):
        viewpoint_node_ids(
            trunk,
            (0.5, 0.5, 0.5),
            np.empty(2, dtype="int32"),
            np.empty(2, dtype="int64"),
        )
