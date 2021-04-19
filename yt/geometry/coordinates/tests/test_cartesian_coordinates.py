# Some tests for the Cartesian coordinates handler

import numpy as np

from yt.testing import assert_equal, fake_amr_ds

# Our canonical tests are that we can access all of our fields and we can
# compute our volume correctly.


def test_cartesian_coordinates():
    # We're going to load up a simple AMR grid and check its volume
    # calculations and path length calculations.
    ds = fake_amr_ds()
    axes = sorted(set(ds.coordinates.axis_name.values()))
    for i, axis in enumerate(axes):
        dd = ds.all_data()
        fi = ("index", axis)
        fd = ("index", f"d{axis}")
        fp = ("index", f"path_element_{axis}")
        ma = np.argmax(dd[fi])
        assert_equal(dd[fi][ma] + dd[fd][ma] / 2.0, ds.domain_right_edge[i])
        mi = np.argmin(dd[fi])
        assert_equal(dd[fi][mi] - dd[fd][mi] / 2.0, ds.domain_left_edge[i])
        assert_equal(dd[fd].min(), ds.index.get_smallest_dx())
        assert_equal(dd[fd].max(), (ds.domain_width / ds.domain_dimensions)[i])
        assert_equal(dd[fd], dd[fp])
    assert_equal(
        dd[("index", "cell_volume")].sum(dtype="float64"), ds.domain_width.prod()
    )
