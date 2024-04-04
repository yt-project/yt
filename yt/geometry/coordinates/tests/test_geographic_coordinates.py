# Some tests for the geographic coordinates handler

import numpy as np
import pytest
from numpy.testing import assert_equal

from yt.testing import assert_rel_equal, fake_amr_ds

# Our canonical tests are that we can access all of our fields and we can
# compute our volume correctly.


@pytest.mark.parametrize("geometry", ("geographic", "internal_geographic"))
def test_geographic_coordinates(geometry):
    # We're going to load up a simple AMR grid and check its volume
    # calculations and path length calculations.

    # Note that we are setting it up to have an altitude of 1000 maximum, which
    # means our volume will be that of a shell 1000 wide, starting at r of
    # whatever our surface_height is set to.
    ds = fake_amr_ds(geometry=geometry)
    if geometry == "geographic":
        ds.surface_height = ds.quan(5000.0, "code_length")
        inner_r = ds.surface_height
        outer_r = ds.surface_height + ds.domain_width[2]
    else:
        ds.outer_radius = ds.quan(5000.0, "code_length")
        inner_r = ds.outer_radius - ds.domain_right_edge[2]
        outer_r = ds.outer_radius
    radial_axis = ds.coordinates.radial_axis
    axes = ds.coordinates.axis_order
    for i, axis in enumerate(axes):
        dd = ds.all_data()
        fi = ("index", axis)
        fd = ("index", f"d{axis}")
        ma = np.argmax(dd[fi])
        assert_equal(dd[fi][ma] + dd[fd][ma] / 2.0, ds.domain_right_edge[i].d)
        mi = np.argmin(dd[fi])
        assert_equal(dd[fi][mi] - dd[fd][mi] / 2.0, ds.domain_left_edge[i].d)
        assert_equal(dd[fd].max(), (ds.domain_width / ds.domain_dimensions)[i].d)

    assert_equal(dd["index", "dtheta"], dd["index", "dlatitude"] * np.pi / 180.0)
    assert_equal(dd["index", "dphi"], dd["index", "dlongitude"] * np.pi / 180.0)
    # Note our terrible agreement here.
    assert_rel_equal(
        dd[("index", "cell_volume")].sum(dtype="float64"),
        (4.0 / 3.0) * np.pi * (outer_r**3 - inner_r**3),
        14,
    )
    assert_equal(
        dd["index", f"path_element_{radial_axis}"], dd["index", f"d{radial_axis}"]
    )
    assert_equal(dd["index", f"path_element_{radial_axis}"], dd["index", "dr"])
    # Note that latitude corresponds to theta, longitude to phi
    assert_equal(
        dd["index", "path_element_latitude"],
        dd["index", "r"] * dd["index", "dlatitude"] * np.pi / 180.0,
    )
    assert_equal(
        dd["index", "path_element_longitude"],
        (
            dd["index", "r"]
            * dd["index", "dlongitude"]
            * np.pi
            / 180.0
            * np.sin((90 - dd["index", "latitude"]) * np.pi / 180.0)
        ),
    )
    # We also want to check that our radius is correct
    offset, factor = ds.coordinates._retrieve_radial_offset()
    radius = factor * dd["index", radial_axis] + offset
    assert_equal(dd["index", "r"], radius)


@pytest.mark.parametrize("geometry", ("geographic", "internal_geographic"))
def test_geographic_conversions(geometry):
    ds = fake_amr_ds(geometry=geometry)
    ad = ds.all_data()
    lats = ad["index", "latitude"]
    dlats = ad["index", "dlatitude"]
    theta = ad["index", "theta"]
    dtheta = ad["index", "dtheta"]

    # check that theta = 0, pi at latitudes of 90, -90
    south_pole_id = np.where(lats == np.min(lats))[0][0]
    north_pole_id = np.where(lats == np.max(lats))[0][0]

    # check that we do in fact have -90, 90 exactly
    assert lats[south_pole_id] - dlats[south_pole_id] / 2.0 == -90.0
    assert lats[north_pole_id] + dlats[north_pole_id] / 2.0 == 90.0

    # check that theta=0 at the north pole, np.pi at the south
    assert theta[north_pole_id] - dtheta[north_pole_id] / 2 == 0.0
    assert theta[south_pole_id] + dtheta[south_pole_id] / 2 == np.pi

    # check that longitude-phi conversions
    phi = ad["index", "phi"]
    dphi = ad["index", "dphi"]
    lons = ad["index", "longitude"]
    dlon = ad["index", "dlongitude"]
    lon_180 = np.where(lons == np.max(lons))[0][0]
    lon_neg180 = np.where(lons == np.min(lons))[0][0]
    # check we have -180, 180 exactly
    assert lons[lon_neg180] - dlon[lon_neg180] / 2.0 == -180.0
    assert lons[lon_180] + dlon[lon_180] / 2.0 == 180.0
    # check that those both convert to phi = np.pi
    assert phi[lon_neg180] - dphi[lon_neg180] / 2.0 == np.pi
    assert phi[lon_180] + dphi[lon_180] / 2.0 == np.pi

    # check that z = +/- radius at +/-90
    # default expected axis order: lat, lon, radial axis
    r_val = ds.coordinates._retrieve_radial_offset()[0]
    coords = np.zeros((2, 3))
    coords[0, 0] = 90.0
    coords[1, 0] = -90.0
    xyz = ds.coordinates.convert_to_cartesian(coords)
    z = xyz[:, 2]
    assert z[0] == r_val
    assert z[1] == -r_val
