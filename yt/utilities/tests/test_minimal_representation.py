import os.path

import yt
from yt.config import ytcfg
from yt.testing import assert_equal, assert_raises, requires_file

G30 = "IsolatedGalaxy/galaxy0030/galaxy0030"


def setup():
    ytcfg["yt", "serialize"] = True


def teardown():
    ytcfg["yt", "serialize"] = False


@requires_file(G30)
def test_store():
    ds = yt.load(G30)
    store = ds.parameter_filename + ".yt"
    field = "density"
    if os.path.isfile(store):
        os.remove(store)

    proj1 = ds.proj(field, "z")
    sp = ds.sphere(ds.domain_center, (4, "kpc"))
    proj2 = ds.proj(field, "z", data_source=sp)

    proj1_c = ds.proj(field, "z")
    assert_equal(proj1[field], proj1_c[field])

    proj2_c = ds.proj(field, "z", data_source=sp)
    assert_equal(proj2[field], proj2_c[field])

    def fail_for_different_method():
        proj2_c = ds.proj(field, "z", data_source=sp, method="max")
        assert_equal(proj2[field], proj2_c[field])

    # A note here: a unyt.exceptions.UnitOperationError is raised
    # and caught by numpy, which reraises a ValueError
    assert_raises(ValueError, fail_for_different_method)

    def fail_for_different_source():
        sp = ds.sphere(ds.domain_center, (2, "kpc"))
        proj2_c = ds.proj(field, "z", data_source=sp, method="integrate")
        assert_equal(proj2_c[field], proj2[field])

    assert_raises(AssertionError, fail_for_different_source)
