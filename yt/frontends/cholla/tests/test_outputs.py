import yt
from yt.frontends.cholla.api import ChollaDataset
from yt.testing import assert_equal, requires_file
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

_fields = (
    ("gas", "temperature"),
    ("gas", "density"),
)

ChollaSimple = "ChollaSimple/0.h5"


@requires_file(ChollaSimple)
def test_ChollaDataset():
    assert isinstance(data_dir_load(ChollaSimple), ChollaDataset)


@requires_file(ChollaSimple)
def test_ChollaSimple_fields():

    expected_fields = [
        "Energy",
        "GasEnergy",
        "density",
        "momentum_x",
        "momentum_y",
        "momentum_z",
        "scalar0",
    ]

    ds = yt.load(ChollaSimple)
    assert_equal(str(ds), "0.h5")
    ad = ds.all_data()

    # Check all the expected fields exist and can be accessed
    for field in expected_fields:
        assert ("cholla", field) in ds.field_list

        # test that field access works
        ad["cholla", field]


@requires_file(ChollaSimple)
def test_ChollaSimple_derived_fields():

    expected_derived_fields = [
        "density",
        "momentum_x",
        "momentum_y",
        "momentum_z",
        "metallicity",
    ]

    ds = yt.load(ChollaSimple)
    ad = ds.all_data()

    # Check all the expected fields exist and can be accessed
    for field in expected_derived_fields:
        assert ("gas", field) in ds.derived_field_list

        # test that field access works
        ad["gas", field]


_fields_chollasimple = (
    ("cholla", "GasEnergy"),
    ("gas", "temperature"),
    ("gas", "density"),
    ("gas", "metallicity"),
)


@requires_ds(ChollaSimple)
def test_cholla_data():
    ds = data_dir_load(ChollaSimple)
    assert_equal(str(ds), "0.h5")
    for test in small_patch_amr(
        ds, _fields_chollasimple, input_center="c", input_weight="ones"
    ):
        test_cholla_data.__name__ = test.description
        yield test
