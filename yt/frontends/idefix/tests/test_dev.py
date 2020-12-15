import os

from yt.frontends.idefix.api import IdefixDataset
from yt.loaders import load
from yt.testing import requires_file

idefix_khi = os.path.join("idefix", "KHI", "dump.0001.dmp")


@requires_file(idefix_khi)
def test_load():
    ds = load(idefix_khi)
    assert isinstance(ds, IdefixDataset)
    assert ds.dimensionality == 2


@requires_file(idefix_khi)
def test_region():
    ds = load(idefix_khi)
    ds.r[:]


@requires_file(idefix_khi)
def test_fields():
    ds = load(idefix_khi)
    expected = [("idefix", "Vc-RHO"), ("idefix", "Vc-VX1"), ("idefix", "Vc-VX2")]
    assert ds.field_list == expected

    assert ("gas", "density") in ds.derived_field_list
    assert ("gas", "velocity_x") in ds.derived_field_list
    assert ("gas", "velocity_y") in ds.derived_field_list


@requires_file(idefix_khi)
def test_get_data():
    ds = load(idefix_khi)
    ds.r[:]["density"]
