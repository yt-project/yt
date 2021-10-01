# Some tests for the Spherical coordinates handler

import numpy as np

from yt import SlicePlot
from yt.testing import add_noise_fields, assert_almost_equal, assert_equal, fake_amr_ds
from yt.utilities.answer_testing.framework import GenericImageTest

# Our canonical tests are that we can access all of our fields and we can
# compute our volume correctly.


def test_spherical_coordinates():
    # We're going to load up a simple AMR grid and check its volume
    # calculations and path length calculations.
    ds = fake_amr_ds(geometry="spherical")
    axes = ["r", "theta", "phi"]
    for i, axis in enumerate(axes):
        dd = ds.all_data()
        fi = ("index", axis)
        fd = ("index", f"d{axis}")
        ma = np.argmax(dd[fi])
        assert_equal(dd[fi][ma] + dd[fd][ma] / 2.0, ds.domain_right_edge[i].d)
        mi = np.argmin(dd[fi])
        assert_equal(dd[fi][mi] - dd[fd][mi] / 2.0, ds.domain_left_edge[i].d)
        assert_equal(dd[fd].max(), (ds.domain_width / ds.domain_dimensions)[i].d)
    # Note that we're using a lot of funny transforms to get to this, so we do
    # not expect to get actual agreement.  This is a bit of a shame, but I
    # don't think it is avoidable as of right now.  Real datasets will almost
    # certainly be correct, if this is correct to 3 decimel places.
    assert_almost_equal(
        dd[("index", "cell_volume")].sum(dtype="float64"),
        (4.0 / 3.0) * np.pi * ds.domain_width[0] ** 3,
    )
    assert_equal(dd["index", "path_element_r"], dd["index", "dr"])
    assert_equal(
        dd["index", "path_element_theta"], dd["index", "r"] * dd["index", "dtheta"]
    )
    assert_equal(
        dd["index", "path_element_phi"],
        (dd["index", "r"] * dd["index", "dphi"] * np.sin(dd["index", "theta"])),
    )


def test_noise_plots():
    ds = fake_amr_ds(geometry="spherical")
    add_noise_fields(ds)

    def create_image(filename_prefix):
        fields = ["noise%d" % i for i in range(4)]

        for normal in ("phi", "theta"):
            p = SlicePlot(ds, normal, fields)
            p.save(f"{filename_prefix}_{normal}")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_noise_plot_lin"
    test_noise_plots.__name__ = test.description
    yield test
