from yt import load_hdf5_file
from yt.testing import assert_almost_equal, assert_equal, requires_file, requires_module

turb_vels = "UnigridData/turb_vels.h5"

_existing_fields = (
    "Bx",
    "By",
    "Bz",
    "Density",
    "MagneticEnergy",
    "Temperature",
    "turb_x-velocity",
    "turb_y-velocity",
    "turb_z-velocity",
    "x-velocity",
    "y-velocity",
    "z-velocity",
)


@requires_file(turb_vels)
@requires_module("h5py")
def test_load_hdf5_file():
    ds1 = load_hdf5_file(turb_vels)
    assert_equal(ds1.domain_dimensions, [256, 256, 256])
    for field_name in _existing_fields:
        assert ("stream", field_name) in ds1.field_list
    assert_equal(ds1.r[:]["ones"].size, 256 * 256 * 256)
    assert_equal(ds1.r[:]["Density"].size, 256 * 256 * 256)
    # Now we test that we get the same results regardless of our decomp
    ds2 = load_hdf5_file(turb_vels, nchunks=19)
    assert_equal(ds2.domain_dimensions, [256, 256, 256])
    assert_equal(ds2.r[:]["ones"].size, 256 * 256 * 256)
    assert_equal(ds2.r[:]["Density"].size, 256 * 256 * 256)
    assert_almost_equal(ds2.r[:]["Density"].min(), ds1.r[:]["Density"].min())
    assert_almost_equal(ds2.r[:]["Density"].max(), ds1.r[:]["Density"].max())
    assert_almost_equal(ds2.r[:]["Density"].std(), ds1.r[:]["Density"].std())
