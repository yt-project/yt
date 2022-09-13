import numpy as np
from nose.tools import assert_equal

from yt.utilities.lib.cykdtree import utils  # type: ignore
from yt.utilities.lib.cykdtree.tests import assert_less_equal, parametrize


def test_max_pts():
    pts = np.arange(5 * 3).reshape((5, 3)).astype("float64")
    out = utils.py_max_pts(pts)
    np.testing.assert_allclose(out, np.max(pts, axis=0))


def test_min_pts():
    pts = np.arange(5 * 3).reshape((5, 3)).astype("float64")
    out = utils.py_min_pts(pts)
    np.testing.assert_allclose(out, np.min(pts, axis=0))


@parametrize(N=(10), ndim=(2, 3), Lidx=(0, 5), Ridx=(5, 9))
def test_argmax_pts_dim(N=10, ndim=2, Lidx=0, Ridx=9):
    d = ndim - 1
    pts = np.random.rand(N, ndim).astype("float64")
    idx = np.argsort(pts[:, d]).astype("uint64")
    out = utils.py_argmax_pts_dim(pts, idx, d, Lidx, Ridx)
    assert_equal(out, np.argmax(pts[idx[Lidx : (Ridx + 1)], d]) + Lidx)


@parametrize(N=(10), ndim=(2, 3), Lidx=(0, 5), Ridx=(5, 9))
def test_argmin_pts_dim(N=10, ndim=2, Lidx=0, Ridx=9):
    d = ndim - 1
    pts = np.random.rand(N, ndim).astype("float64")
    idx = np.argsort(pts[:, d]).astype("uint64")
    out = utils.py_argmin_pts_dim(pts, idx, d, Lidx, Ridx)
    assert_equal(out, np.argmin(pts[idx[Lidx : (Ridx + 1)], d]) + Lidx)


@parametrize(N=(0, 10, 11), ndim=(2, 3))
def test_quickSort(N=10, ndim=2):
    d = ndim - 1
    np.random.seed(10)
    pts = np.random.rand(N, ndim).astype("float64")
    idx = utils.py_quickSort(pts, d)
    assert_equal(idx.size, N)
    if N != 0:
        np.testing.assert_allclose(idx, np.argsort(pts[:, d]))


@parametrize(N=(0, 10, 11), ndim=(2, 3))
def test_insertSort(N=10, ndim=2):
    d = ndim - 1
    np.random.seed(10)
    pts = np.random.rand(N, ndim).astype("float64")
    idx = utils.py_insertSort(pts, d)
    assert_equal(idx.size, N)
    if N != 0:
        np.testing.assert_allclose(idx, np.argsort(pts[:, d]))


@parametrize(N=(0, 10, 11), ndim=(2, 3))
def test_pivot(N=10, ndim=2):
    d = ndim - 1
    np.random.seed(10)
    pts = np.random.rand(N, ndim).astype("float64")
    q, idx = utils.py_pivot(pts, d)
    if N == 0:
        np.testing.assert_equal(q, -1)
    else:
        piv = pts[idx[q], d]
        nmax = 7 * N / 10 + 6
        assert_less_equal(np.sum(pts[:, d] < piv), nmax)
        assert_less_equal(np.sum(pts[:, d] > piv), nmax)


@parametrize(N=(0, 10, 11), ndim=(2, 3))
def test_partition(N=10, ndim=2):
    d = ndim - 1
    p = 0
    np.random.seed(10)
    pts = np.random.rand(N, ndim).astype("float64")
    q, idx = utils.py_partition(pts, d, p)
    if N == 0:
        assert_equal(q, -1)
    else:
        piv = pts[p, d]
        np.testing.assert_approx_equal(pts[idx[q], d], piv)
        np.testing.assert_array_less(pts[idx[:q], d], piv)
        np.testing.assert_array_less(piv, pts[idx[(q + 1) :], d])


@parametrize(N=(0, 10, 11), ndim=(2, 3))
def test_partition_given_pivot(N=10, ndim=2):
    d = ndim - 1
    np.random.seed(10)
    pts = np.random.rand(N, ndim).astype("float64")
    if N == 0:
        piv_list = [0.5]
    else:
        piv_list = [0.5, np.median(pts[:, d])]
    for piv in piv_list:
        q, idx = utils.py_partition_given_pivot(pts, d, piv)
        if N == 0:
            assert_equal(q, -1)
        else:
            assert_less_equal(pts[idx[q], d], piv)
            np.testing.assert_array_less(pts[idx[:q], d], piv)
            np.testing.assert_array_less(piv, pts[idx[(q + 1) :], d])


@parametrize(N=(0, 10, 11), ndim=(2, 3))
def test_select(N=10, ndim=2):
    d = ndim - 1
    np.random.seed(10)
    pts = np.random.rand(N, ndim).astype("float64")
    p = int(N) // 2 + int(N) % 2
    q, idx = utils.py_select(pts, d, p)
    assert_equal(idx.size, N)
    if N == 0:
        assert_equal(q, -1)
    else:
        assert_equal(q, p - 1)
        med = np.median(pts[:, d])
        np.testing.assert_array_less(pts[idx[:q], d], med)
        np.testing.assert_array_less(med, pts[idx[(q + 1) :], d])
        if N % 2:
            np.testing.assert_approx_equal(pts[idx[q], d], med)
        else:
            np.testing.assert_array_less(pts[idx[q], d], med)


@parametrize(N=(0, 10, 11), ndim=(2, 3), use_sliding_midpoint=(False, True))
def test_split(N=10, ndim=2, use_sliding_midpoint=False):
    np.random.seed(10)
    pts = np.random.rand(N, ndim).astype("float64")
    p = int(N) // 2 + int(N) % 2
    q, d, idx = utils.py_split(pts, use_sliding_midpoint=use_sliding_midpoint)
    assert_equal(idx.size, N)
    if N == 0:
        assert_equal(q, -1)
    else:
        if use_sliding_midpoint:
            # Midpoint
            med = 0.5 * (np.min(pts[:, d]) + np.max(pts[:, d]))
            np.testing.assert_array_less(pts[idx[:q], d], med)
            np.testing.assert_array_less(med, pts[idx[(q + 1) :], d])
            np.testing.assert_array_less(pts[idx[q], d], med)
            # Sliding midpoint (slide to minimum)
            q, d, idx = utils.py_split(
                pts,
                mins=-1 * np.ones(ndim),
                maxs=np.ones(ndim),
                use_sliding_midpoint=True,
            )
            med = np.min(pts[:, d])
            assert_equal(q, 0)
            np.testing.assert_array_less(pts[idx[:q], d], med)
            np.testing.assert_array_less(med, pts[idx[(q + 1) :], d])
            np.testing.assert_approx_equal(pts[idx[q], d], med)
            # Sliding midpoint (slide to maximum)
            q, d, idx = utils.py_split(
                pts,
                mins=np.zeros(ndim),
                maxs=2 * np.ones(ndim),
                use_sliding_midpoint=True,
            )
            med = np.max(pts[:, d])
            assert_equal(q, N - 2)
            np.testing.assert_array_less(pts[idx[: (q + 1)], d], med)
            np.testing.assert_approx_equal(pts[idx[q + 1], d], med)
        else:
            assert_equal(q, p - 1)
            med = np.median(pts[:, d])
            np.testing.assert_array_less(pts[idx[:q], d], med)
            np.testing.assert_array_less(med, pts[idx[(q + 1) :], d])
            if N % 2:
                np.testing.assert_approx_equal(pts[idx[q], d], med)
            else:
                np.testing.assert_array_less(pts[idx[q], d], med)
