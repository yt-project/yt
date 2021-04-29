from yt.frontends.gdf.api import GDFDataset
from yt.testing import assert_equal, requires_file, units_override_check
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

_fields = [("gas", "density"), ("gas", "velocity_x")]

sedov = "sedov/sedov_tst_0004.h5"


@requires_ds(sedov)
def test_sedov_tunnel():
    ds = data_dir_load(sedov)
    assert_equal(str(ds), "sedov_tst_0004")
    for test in small_patch_amr(ds, _fields):
        test_sedov_tunnel.__name__ = test.description
        yield test


@requires_file(sedov)
def test_GDFDataset():
    assert isinstance(data_dir_load(sedov), GDFDataset)


@requires_file(sedov)
def test_units_override():
    units_override_check(sedov)
