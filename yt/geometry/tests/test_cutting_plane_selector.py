import numpy as np
import pytest

from yt.geometry.selection_routines import cutting_mixed_spherical_selector


class HelpfulPlaneObject:
    # a bare-bones skeleton of a data object to use for initializing
    # the cython cutting plane selector
    def __init__(self, normal, plane_center):
        self._d = -1 * np.dot(normal, plane_center)
        self._norm_vec = normal


@pytest.fixture
def xy_plane_at_001():
    normal = np.array([0.0, 0.0, 1.0])
    plane_center = np.array([0.0, 0.0, 1.0])
    return HelpfulPlaneObject(normal, plane_center)


def test_spherical_cutting_plane_spots(xy_plane_at_001):
    # a couple of manual spot checks for intersection of a plane
    # with some spherical volume elements

    # initialize selector
    scp = cutting_mixed_spherical_selector(xy_plane_at_001)
    assert scp.r_min == 1.0

    # left/right edge values are given in spherical coordinates with
    # order of (r, theta, phi) where
    #   theta is the azimuthal/latitudinal
    #   phi is the polar/longitudinal angle (bounds 0 to 2pi).

    def _in_rads(x):
        return x * np.pi / 180

    # should intersect
    left_edge = np.array([0.8, _in_rads(5), _in_rads(5)])
    right_edge = np.array([1.2, _in_rads(45), _in_rads(45)])
    assert scp._select_single_bbox(left_edge, right_edge)

    # should not intersect
    left_edge = np.array([0.1, _in_rads(90), _in_rads(5)])
    right_edge = np.array([0.4, _in_rads(110), _in_rads(45)])
    assert scp._select_single_bbox(left_edge, right_edge) == 0


def test_spherical_cutting_plane(xy_plane_at_001):
    import numpy as np

    from yt.testing import fake_amr_ds

    ds = fake_amr_ds(geometry="spherical")

    # this plane will miss the dataset entirely
    normal = np.array([0.0, 0.0, 1.0])
    plane_center = np.array([0.0, 0.0, 1.1])
    slc = ds.cutting_mixed(normal, plane_center)
    assert len(slc[("stream", "Density")]) == 0

    # this one will not.
    normal = np.array([0.0, 0.0, 1.0])
    plane_center = np.array([0.0, 0.0, 0.5])
    slc = ds.cutting_mixed(normal, plane_center)
    r = slc[("index", "r")]
    theta = slc[("index", "theta")]
    # r cannot be smaller than the distance from plane to origin
    assert np.min(r) >= plane_center[2]

    # how close the z value is to the plane's z value
    # depends on the size of the elements, the closeness here
    # was found empirically
    z = r * np.cos(theta)
    max_z = np.max(np.abs(z.to("code_length").d - 0.5))
    assert np.isclose(max_z, 0.04212724)


def test_spherical_cutting_plane_frb():
    # future image test
    pass
