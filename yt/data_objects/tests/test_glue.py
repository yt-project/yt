import nose

from yt.config import ytcfg
from yt.testing import fake_random_ds, requires_module


def setup_func():
    ytcfg["yt", "internals", "within_testing"] = True


@requires_module("glue.core")
@nose.with_setup(setup_func)
def test_glue_data_object():
    ds = fake_random_ds(16)
    ad = ds.all_data()
    ad.to_glue([("gas", "density")])
