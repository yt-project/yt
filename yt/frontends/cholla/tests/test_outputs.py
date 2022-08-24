import yt
from yt.frontends.ramses.api import ChollaDataset
from yt.testing import assert_equal, requires_file
from yt.utilities.answer_testing.framework import data_dir_load, requires_ds

_fields = (
    ("gas", "temperature"),
    ("gas", "density"),
)

cholla_0 = "cholla_0.h5"


@requires_ds(cholla_0)
def test_cholla_0():

    expected_fields = [
        "Density",
        "x-velocity",
        "y-velocity",
        "z-velocity",
        "Pres_IR",
        "Pressure",
        "Metallicity",
        "HII",
        "HeII",
        "HeIII",
    ]

    ds = yt.load(cholla_0)
    assert_equal(str(ds), "cholla_0")
    yt.ProjectionPlot(ds, "x", ("gas", "density"))
    ad = ds.all_data()

    # Check all the expected fields exist and can be accessed
    for f in expected_fields:
        assert f in ds.derived_field_list
        ad[f]

    for field in expected_fields:
        assert ("cholla", field) in ds.field_list

        # test that field access works
        ad["cholla", field]


@requires_file(cholla_0)
def test_ChollaDataset():
    assert isinstance(data_dir_load(cholla_0), ChollaDataset)
