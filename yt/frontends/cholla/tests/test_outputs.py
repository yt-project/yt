import yt
from yt.frontends.cholla.api import ChollaDataset
from yt.testing import assert_equal, requires_file
from yt.utilities.answer_testing.framework import data_dir_load

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
