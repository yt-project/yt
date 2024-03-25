import numpy as np
import pytest

from yt import load_uniform_grid
from yt.geometry.selection_routines import cutting_mixed_spherical_selector
from yt.testing import fake_amr_ds


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


def _in_rads(x):
    return x * np.pi / 180


def test_spherical_cutting_plane_spots(xy_plane_at_001):
    # a couple of manual spot checks for intersection of a plane
    # with some spherical volume elements

    # initialize selector
    scp = cutting_mixed_spherical_selector(xy_plane_at_001)
    assert scp.r_min == 1.0

    # left/right edge values are given in spherical coordinates with
    # order of (r, theta, phi) where
    #   theta is the colatitude (bounds 0 to pi)
    #   phi is the azimuth (bounds 0 to 2pi).

    # should intersect
    left_edge = np.array([0.8, _in_rads(5), _in_rads(5)])
    right_edge = np.array([1.2, _in_rads(45), _in_rads(45)])
    assert scp._select_single_bbox(left_edge, right_edge)

    # should not intersect
    left_edge = np.array([0.1, _in_rads(90), _in_rads(5)])
    right_edge = np.array([0.4, _in_rads(110), _in_rads(45)])
    assert scp._select_single_bbox(left_edge, right_edge) == 0


def test_large_angular_range():
    # check that large elements are still selected

    # these edges define a single element that is a spherical shell of finite
    # thickness spanning a hemisphere. The bounds of the element all fall on
    # one side of the test plane, so these checks rely on the additional angular
    # verts that get added for large elements
    left_edge = np.array([0.8, 0.01, 0.01])
    right_edge = np.array([1.0, np.pi - 0.01, np.pi - 0.01])

    for y_pos in np.linspace(0.1, 0.99, 10):
        normal = np.array([0.0, 1.0, 0.0])
        plane_center = np.array([0.0, y_pos, 0.0])
        xz_plane = HelpfulPlaneObject(normal, plane_center)
        scp = cutting_mixed_spherical_selector(xz_plane)

        selected = scp._select_single_bbox(left_edge, right_edge)
        assert selected

        lev = np.array(
            [
                [
                    0,
                ]
            ],
            dtype=np.int32,
        )
        left_edges = np.array(
            [
                left_edge,
            ]
        )
        right_edges = np.array(
            [
                right_edge,
            ]
        )
        grid_sel = scp.select_grids(left_edges, right_edges, lev)
        assert grid_sel


def test_large_angular_range_ds():
    # checks that a ds in spherical coords with a single grid spanning
    # a large angular range gets selected properly.
    bbox = np.array([[0.5, 1.0], [0, np.pi], [0, np.pi]])

    shp = (32,) * 3
    data = {"density": np.random.random(shp)}

    ds = load_uniform_grid(
        data,
        shp,
        bbox=bbox,
        geometry="spherical",
        axis_order=("r", "theta", "phi"),
        length_unit="m",
        nprocs=1,
    )

    normal = ds.arr([0.0, 1.0, 0], "code_length")
    center = ds.arr([0.0, 0.2, 0.0], "code_length")
    slc = ds.cutting_mixed(normal, center)

    le = ds.index.grid_left_edge
    re = ds.index.grid_right_edge
    lev = ds.index.grid_levels
    selected = slc.selector.select_grids(le, re, lev)
    assert np.all(selected)


def test_spherical_cutting_plane(xy_plane_at_001):

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
