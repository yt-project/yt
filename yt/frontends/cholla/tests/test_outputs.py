import yt
from yt.frontends.cholla.api import ChollaDataset
from yt.testing import assert_equal, requires_file
from yt.utilities.answer_testing.framework import data_dir_load

_fields = (
    ("gas", "temperature"),
    ("gas", "density"),
)

cholla_0 = "cholla_0.h5"


@requires_file(cholla_0)
def test_ChollaDataset():
    assert isinstance(data_dir_load(cholla_0), ChollaDataset)


@requires_file(cholla_0)
def test_cholla_0_fields():

    expected_fields = [
        "Energy",
        "GasEnergy",
        "density",
        "momentum_x",
        "momentum_y",
        "momentum_z",
        "scalar0",
    ]

    ds = yt.load(cholla_0)
    assert_equal(str(ds), "cholla_0.h5")
    ad = ds.all_data()

    # Check all the expected fields exist and can be accessed
    for field in expected_fields:
        assert ("cholla", field) in ds.field_list
        # assert ("gas", field) in ds.derived_field_list

        # test that field access works
        ad["cholla", field]


@requires_file(cholla_0)
def test_cholla_0_derived_fields():

    expected_derived_fields = [
        "density",
        "momentum_x",
        "momentum_y",
        "momentum_z",
        "metallicity",
    ]

    ds = yt.load(cholla_0)
    ad = ds.all_data()

    # Check all the expected fields exist and can be accessed
    for field in expected_derived_fields:
        assert ("gas", field) in ds.derived_field_list

        # test that field access works
        ad["gas", field]
