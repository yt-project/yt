from yt.frontends.chombo.api import ChomboDataset, Orion2Dataset, PlutoDataset
from yt.testing import assert_equal, requires_file, units_override_check
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

_fields = (
    ("gas", "density"),
    ("gas", "velocity_magnitude"),
    ("gas", "magnetic_field_x"),
)

gc = "GaussianCloud/data.0077.3d.hdf5"


@requires_ds(gc)
def test_gc():
    ds = data_dir_load(gc)
    assert_equal(str(ds), "data.0077.3d.hdf5")
    for test in small_patch_amr(ds, _fields):
        test_gc.__name__ = test.description
        yield test


tb = "TurbBoxLowRes/data.0005.3d.hdf5"


@requires_ds(tb)
def test_tb():
    ds = data_dir_load(tb)
    assert_equal(str(ds), "data.0005.3d.hdf5")
    for test in small_patch_amr(ds, _fields):
        test_tb.__name__ = test.description
        yield test


iso = "IsothermalSphere/data.0000.3d.hdf5"


@requires_ds(iso)
def test_iso():
    ds = data_dir_load(iso)
    assert_equal(str(ds), "data.0000.3d.hdf5")
    for test in small_patch_amr(ds, _fields):
        test_iso.__name__ = test.description
        yield test


_zp_fields = (("chombo", "rhs"), ("chombo", "phi"))
zp = "ZeldovichPancake/plt32.2d.hdf5"


@requires_ds(zp)
def test_zp():
    ds = data_dir_load(zp)
    assert_equal(str(ds), "plt32.2d.hdf5")
    for test in small_patch_amr(ds, _zp_fields, input_center="c", input_weight="rhs"):
        test_zp.__name__ = test.description
        yield test


kho = "KelvinHelmholtz/data.0004.hdf5"


@requires_ds(kho)
def test_kho():
    ds = data_dir_load(kho)
    assert_equal(str(ds), "data.0004.hdf5")
    for test in small_patch_amr(ds, _fields):
        test_kho.__name__ = test.description
        yield test


@requires_file(zp)
def test_ChomboDataset():
    assert isinstance(data_dir_load(zp), ChomboDataset)


@requires_file(gc)
def test_Orion2Dataset():
    assert isinstance(data_dir_load(gc), Orion2Dataset)


@requires_file(kho)
def test_PlutoDataset():
    assert isinstance(data_dir_load(kho), PlutoDataset)


@requires_file(zp)
def test_units_override_zp():
    units_override_check(zp)


@requires_file(gc)
def test_units_override_gc():
    units_override_check(gc)


@requires_file(kho)
def test_units_override_kho():
    units_override_check(kho)
