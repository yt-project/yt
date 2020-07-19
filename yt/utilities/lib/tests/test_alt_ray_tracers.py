"""Tests for non-cartesian ray tracers."""
import numpy as np

from yt.testing import amrspace, assert_equal, assert_less_equal, assert_true
from yt.utilities.lib.alt_ray_tracers import _cyl2cart, cylindrical_ray_trace

left_grid = right_grid = amr_levels = center_grid = data = None


def setup():
    # set up some sample cylindrical grid data, radiating out from center
    global left_grid, right_grid, amr_levels, center_grid, data
    np.seterr(all="ignore")
    l1, r1, lvl1 = amrspace([0.0, 1.0, 0.0, -1.0, 0.0, 2 * np.pi], levels=(7, 7, 0))
    l2, r2, lvl2 = amrspace([0.0, 1.0, 0.0, 1.0, 0.0, 2 * np.pi], levels=(7, 7, 0))
    left_grid = np.concatenate([l1, l2], axis=0)
    right_grid = np.concatenate([r1, r2], axis=0)
    amr_levels = np.concatenate([lvl1, lvl2], axis=0)
    center_grid = (left_grid + right_grid) / 2.0
    data = np.cos(np.sqrt(np.sum(center_grid[:, :2] ** 2, axis=1))) ** 2  # cos^2


point_pairs = np.array(
    [
        # p1               p2
        ([0.5, -1.0, 0.0], [1.0, 1.0, 0.75 * np.pi]),  # Everything different
        ([0.5, -1.0, 0.0], [0.5, 1.0, 0.75 * np.pi]),  # r same
        ([0.5, -1.0, 0.0], [0.5, 1.0, np.pi]),  # diagonal through z-axis
        # straight through z-axis
        ([0.5, 0.0, 0.0], [0.5, 0.0, np.pi]),
        # ([0.5, 0.0, np.pi*3/2 + 0.0], [0.5, 0.0, np.pi*3/2 + np.pi]),
        # ([0.5, 0.0, np.pi/2 + 0.0], [0.5, 0.0, np.pi/2 + np.pi]),
        # ([0.5, 0.0, np.pi + 0.0], [0.5, 0.0, np.pi + np.pi]),
        # const z, not through z-axis
        ([0.5, 0.1, 0.0], [0.5, 0.1, 0.75 * np.pi]),
        # ([0.5, 0.1, np.pi + 0.0], [0.5, 0.1, np.pi + 0.75*np.pi]),
        # ([0.5, 0.1, np.pi*3/2 + 0.0], [0.5, 0.1, np.pi*3/2 + 0.75*np.pi]),
        # ([0.5, 0.1, np.pi/2 + 0.0], [0.5, 0.1, np.pi/2 + 0.75*np.pi]),
        # ([0.5, 0.1, 2*np.pi + 0.0], [0.5, 0.1, 2*np.pi + 0.75*np.pi]),
        # ([0.5, 0.1, np.pi/4 + 0.0], [0.5, 0.1, np.pi/4 + 0.75*np.pi]),
        # ([0.5, 0.1, np.pi*3/8 + 0.0], [0.5, 0.1, np.pi*3/8 + 0.75*np.pi]),
        (
            [0.5, -1.0, 0.75 * np.pi],
            [1.0, 1.0, 0.75 * np.pi],
        ),  # r,z different - theta same
        ([0.5, -1.0, 0.75 * np.pi], [0.5, 1.0, 0.75 * np.pi]),  # z-axis parallel
        ([0.0, -1.0, 0.0], [0.0, 1.0, 0.0]),  # z-axis itself
    ]
)


def check_monotonic_inc(arr):
    assert_true(np.all(0.0 <= (arr[1:] - arr[:-1])))


def check_bounds(arr, blower, bupper):
    assert_true(np.all(blower <= arr))
    assert_true(np.all(bupper >= arr))


def test_cylindrical_ray_trace():
    for pair in point_pairs:
        p1, p2 = pair
        p1cart, p2cart = _cyl2cart(pair)
        pathlen = np.sqrt(np.sum((p2cart - p1cart) ** 2))

        t, s, rztheta, inds = cylindrical_ray_trace(p1, p2, left_grid, right_grid)
        npoints = len(t)

        check_monotonic_inc(t)
        assert_less_equal(0.0, t[0])
        assert_less_equal(t[-1], 1.0)

        check_monotonic_inc(s)
        assert_less_equal(0.0, s[0])
        assert_less_equal(s[-1], pathlen)
        assert_equal(npoints, len(s))

        assert_equal((npoints, 3), rztheta.shape)
        check_bounds(rztheta[:, 0], 0.0, 1.0)
        check_bounds(rztheta[:, 1], -1.0, 1.0)
        check_bounds(rztheta[:, 2], 0.0, 2 * np.pi)
        check_monotonic_inc(rztheta[:, 2])

        assert_equal(npoints, len(inds))
        check_bounds(inds, 0, len(left_grid) - 1)
