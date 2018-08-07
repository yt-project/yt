import numpy as np

from yt.testing import \
    fake_amr_ds, \
    assert_equal, \
    assert_rel_equal

def test_alt_az_coordinates():
    # We're going to load up a simple AMR grid and check its volume
    # calculations and path length calculations.

    # Note that we are setting it up to have an altitude of 1000 maximum, which
    # means our volume will be that of a shell 1000 wide, starting at r of
    # whatever our surface_height is set to.
    ds = fake_amr_ds(geometry="alt_az")
    axes = ["range", "azimuth", "elevation"]
    for i, axis in enumerate(axes):
        dd = ds.all_data()
        fi = ("index", axis)
        fd = ("index", "d%s" % axis)
        ma = np.argmax(dd[fi])
        assert_equal(dd[fi][ma] + dd[fd][ma] / 2.0, ds.domain_right_edge[i].d)
        mi = np.argmin(dd[fi])
        assert_equal(dd[fi][mi] - dd[fd][mi] / 2.0, ds.domain_left_edge[i].d)
        assert_equal(dd[fd].max(), (ds.domain_width/ds.domain_dimensions)[i].d)
    outer_r = ds.domain_width[0]
    assert_equal(dd["index","dtheta"], dd["index","delevation"]*np.pi/180.0)
    assert_equal(dd["index","dphi"], dd["index","dazimuth"]*np.pi/180.0)
    assert_rel_equal(dd["cell_volume"].sum(dtype="float64"),
                     (4.0/3.0) * np.pi * outer_r**3, 10)
    assert_equal(dd["index", "path_element_range"], dd["index", "drange"])
    assert_equal(dd["index", "path_element_range"], dd["index", "dr"])
    # Note that elevation corresponds to theta, azimuth to phi
    assert_equal(dd["index", "path_element_elevation"],
                 dd["index", "r"] * dd["index", "delevation"] * np.pi/180.0)
    assert_equal(dd["index", "path_element_azimuth"],
                 (dd["index", "r"] * dd["index", "dazimuth"] * np.pi/180.0 *
                  np.sin((dd["index", "elevation"] + 90.0) * np.pi/180.0)))
    # We also want to check that our radius is correct
    assert_equal(dd["index","r"], dd["index","range"])
