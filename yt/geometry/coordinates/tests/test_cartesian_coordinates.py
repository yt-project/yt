# Some tests for the Cartesian coordinates handler

import numpy as np

from yt.testing import \
    fake_amr_ds, \
    assert_equal

# Our canonical tests are that we can access all of our fields and we can
# compute our volume correctly.

def test_cartesian_coordinates():
    # We're going to load up a simple AMR grid and check its volume
    # calculations and path length calculations.
    ds = fake_amr_ds()
    axes = list(set(ds.coordinates.axis_name.values()))
    axes.sort()
    for i, axis in enumerate(axes):
        dd = ds.all_data()
        fi = ("index", axis)
        fd = ("index", "d%s" % axis)
        fp = ("index", "path_element_%s" % axis)
        ma = np.argmax(dd[fi])
        yield assert_equal, dd[fi][ma] + dd[fd][ma] / 2.0, ds.domain_right_edge[i]
        mi = np.argmin(dd[fi])
        yield assert_equal, dd[fi][mi] - dd[fd][mi] / 2.0, ds.domain_left_edge[i]
        yield assert_equal, dd[fd].min(), ds.index.get_smallest_dx()
        yield assert_equal, dd[fd].max(), (ds.domain_width/ds.domain_dimensions)[i]
        yield assert_equal, dd[fd], dd[fp]
    yield assert_equal, dd["cell_volume"].sum(dtype="float64"), \
                        ds.domain_width.prod()
