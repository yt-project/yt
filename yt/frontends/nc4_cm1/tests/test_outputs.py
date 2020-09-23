from yt.frontends.nc4_cm1.api import CM1Dataset
from yt.testing import assert_equal, requires_file, units_override_check
from yt.utilities.answer_testing.framework import (
    FieldValuesTest,
    GridValuesTest,
    can_run_ds,
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

_fields = ("thrhopert", "zvort")
cm1sim = "testyt.05500.000000.nc"


@requires_ds(cm1sim, big_data=True)
def test_mesh():
    ds = data_dir_load(cm1sim)
    assert_equal(str(ds), "testyt.05500.000000.nc")

    # run the small_patch_amr tests on safe fields
    ic = ds.domain_center
    for test in small_patch_amr(ds, _fields, input_center=ic, input_weight=None):
        test_mesh.__name__ = test.description
        yield test

    # manually run the Grid and Field Values tests on dbz (do not want to run the
    # ProjectionValuesTest for this field)
    if can_run_ds(ds):
        dso = [None, ("sphere", (ic, (0.1, "unitary")))]
        for field in ["dbz"]:
            yield GridValuesTest(ds, field)
            for axis in [0, 1, 2]:
                for dobj_name in dso:
                    yield FieldValuesTest(ds, field, dobj_name)


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
