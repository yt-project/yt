import numpy as np

from yt.testing import assert_almost_equal, fake_random_ds
from yt.utilities.math_utils import euclidean_dist, periodic_dist


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def test_periodicity():
    # First test the simple case were we find the distance between two points
    a = [0.1, 0.1, 0.1]
    b = [0.9, 0.9, 0.9]
    period = 1.0
    dist = periodic_dist(a, b, period)
    assert_almost_equal(dist, 0.34641016151377535)
    dist = periodic_dist(a, b, period, (True, False, False))
    assert_almost_equal(dist, 1.1489125293076059)
    dist = periodic_dist(a, b, period, (False, True, False))
    assert_almost_equal(dist, 1.1489125293076059)
    dist = periodic_dist(a, b, period, (False, False, True))
    assert_almost_equal(dist, 1.1489125293076059)
    dist = periodic_dist(a, b, period, (True, True, False))
    assert_almost_equal(dist, 0.84852813742385713)
    dist = periodic_dist(a, b, period, (True, False, True))
    assert_almost_equal(dist, 0.84852813742385713)
    dist = periodic_dist(a, b, period, (False, True, True))
    assert_almost_equal(dist, 0.84852813742385713)
    dist = euclidean_dist(a, b)
    assert_almost_equal(dist, 1.3856406460551021)

    # Now test the more complicated cases where we're calculating radii based
    # on data objects
    ds = fake_random_ds(64)

    # First we test flattened data
    data = ds.all_data()
    positions = np.array([data[("index", ax)] for ax in "xyz"])
    c = [0.1, 0.1, 0.1]
    n_tup = tuple(1 for i in range(positions.ndim - 1))
    center = np.tile(
        np.reshape(np.array(c), (positions.shape[0],) + n_tup),
        (1,) + positions.shape[1:],
    )

    dist = periodic_dist(positions, center, period, ds.periodicity)
    assert_almost_equal(dist.min(), 0.00270632938683)
    assert_almost_equal(dist.max(), 0.863319074398)

    dist = euclidean_dist(positions, center)
    assert_almost_equal(dist.min(), 0.00270632938683)
    assert_almost_equal(dist.max(), 1.54531407988)

    # Then grid-like data
    data = ds.index.grids[0]
    positions = np.array([data[("index", ax)] for ax in "xyz"])
    c = [0.1, 0.1, 0.1]
    n_tup = tuple(1 for i in range(positions.ndim - 1))
    center = np.tile(
        np.reshape(np.array(c), (positions.shape[0],) + n_tup),
        (1,) + positions.shape[1:],
    )

    dist = periodic_dist(positions, center, period, ds.periodicity)
    assert_almost_equal(dist.min(), 0.00270632938683)
    assert_almost_equal(dist.max(), 0.863319074398)

    dist = euclidean_dist(positions, center)
    assert_almost_equal(dist.min(), 0.00270632938683)
    assert_almost_equal(dist.max(), 1.54531407988)
