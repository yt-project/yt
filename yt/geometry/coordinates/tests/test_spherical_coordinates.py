# Some tests for the Cylindrical coordinates handler

import numpy as np

from yt.testing import \
    fake_amr_ds, \
    assert_equal, \
    assert_almost_equal

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
        fd = ("index", "d%s" % axis)
        ma = np.argmax(dd[fi])
        yield assert_equal, dd[fi][ma] + dd[fd][ma] / 2.0, ds.domain_right_edge[i].d
        mi = np.argmin(dd[fi])
        yield assert_equal, dd[fi][mi] - dd[fd][mi] / 2.0, ds.domain_left_edge[i].d
        yield assert_equal, dd[fd].max(), (ds.domain_width/ds.domain_dimensions)[i].d
    # Note that we're using a lot of funny transforms to get to this, so we do
    # not expect to get actual agreement.  This is a bit of a shame, but I
    # don't think it is avoidable as of right now.  Real datasets will almost
    # certainly be correct, if this is correct to 3 decimel places.
    yield assert_almost_equal, dd["cell_volume"].sum(dtype="float64"), \
                        (4.0/3.0) * np.pi*ds.domain_width[0]**3, \
                        3 # Allow for GENEROUS error on the volume ...
    yield assert_equal, dd["index", "path_element_r"], dd["index", "dr"]
    yield assert_equal, dd["index", "path_element_theta"], \
                        dd["index", "r"] * dd["index", "dtheta"]
    yield assert_equal, dd["index", "path_element_phi"], \
                        dd["index", "r"] * dd["index", "dphi"] * \
                          np.sin(dd["index","theta"])
