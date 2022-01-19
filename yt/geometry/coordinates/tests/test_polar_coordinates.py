# Some tests for the polar coordinates handler
# (Pretty similar to cylindrical, but different ordering)

import numpy as np

from yt.testing import assert_almost_equal, assert_equal, fake_amr_ds

# Our canonical tests are that we can access all of our fields and we can
# compute our volume correctly.


def test_cylindrical_coordinates():
    # We're going to load up a simple AMR grid and check its volume
    # calculations and path length calculations.
    ds = fake_amr_ds(geometry="polar")
    axes = ["r", "theta", "z"]
    for i, axis in enumerate(axes):
        dd = ds.all_data()
        fi = ("index", axis)
        fd = ("index", f"d{axis}")
        ma = np.argmax(dd[fi])
        assert_equal(dd[fi][ma] + dd[fd][ma] / 2.0, ds.domain_right_edge[i].d)
        mi = np.argmin(dd[fi])
        assert_equal(dd[fi][mi] - dd[fd][mi] / 2.0, ds.domain_left_edge[i].d)
        assert_equal(dd[fd].max(), (ds.domain_width / ds.domain_dimensions)[i].d)
    assert_almost_equal(
        dd[("index", "cell_volume")].sum(dtype="float64"),
        np.pi * ds.domain_width[0] ** 2 * ds.domain_width[2],
    )
    assert_equal(dd["index", "path_element_r"], dd["index", "dr"])
    assert_equal(dd["index", "path_element_z"], dd["index", "dz"])
    assert_equal(
        dd["index", "path_element_theta"], dd["index", "r"] * dd["index", "dtheta"]
    )
