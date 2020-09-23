from yt.frontends.nc4_cm1.api import CM1Dataset
from yt.testing import assert_equal, requires_file, units_override_check
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

_fields = ("dbz", "thrhopert", "zvort")
cm1sim = "testyt.05500.000000.nc"


@requires_ds(cm1sim, big_data=True)
def test_mesh():
    ds = data_dir_load(cm1sim)
    print(ds)
    assert_equal(str(ds), "testyt.05500.000000.nc")
    for test in small_patch_amr(ds, _fields):
        test_mesh.__name__ = test.description
        yield test


@requires_file(cm1sim)
def test_CM1Dataset():
    assert isinstance(data_dir_load(cm1sim), CM1Dataset)


@requires_file(cm1sim)
def test_units_override():
    units_override_check(cm1sim)


@requires_file(cm1sim)
def test_dims_and_meta():
    ds = data_dir_load(cm1sim)

    known_dims = ["time", "zf", "zh", "yf", "yh", "xf", "xh"]
    dims = ds.parameters["coords"].dims.keys()

    ## check the file for 2 grids and a time dimension -
    ## (time, xf, xh, yf, yh, zf, zh). The dimesions ending in
    ## f are the staggered velocity grid components and the
    ## dimensions ending in h are the scalar grid components
    assert_equal(len(dims), 7)
    for kdim in known_dims:
        assert kdim in dims

    ## check the simulation time
    assert_equal(ds.parameters["time"], 5500.0)


@requires_file(cm1sim)
def test_if_cm1():
    ds = data_dir_load(cm1sim)
    assert float(ds.parameters["lofs_version"]) >= 1.0
