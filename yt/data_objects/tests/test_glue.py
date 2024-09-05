import pytest

from yt.config import ytcfg
from yt.testing import fake_random_ds, requires_module_pytest as requires_module


@pytest.fixture
def within_testing():
    old_value = ytcfg["yt", "internals", "within_testing"]
    ytcfg["yt", "internals", "within_testing"] = True
    yield
    ytcfg["yt", "internals", "within_testing"] = old_value


@pytest.mark.usefixtures("within_testing")
@requires_module("glue", "astropy")
def test_glue_data_object():
    ds = fake_random_ds(16)
    ad = ds.all_data()
    ad.to_glue([("gas", "density")])
