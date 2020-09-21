from collections import OrderedDict

import os
import stat
import weakref
import numpy as np

from yt.utilities.on_demand_imports import _xarray as xarray
from yt.frontends.nc4_cm1.api import CM1Dataset
from yt.testing import (
    assert_equal,
    requires_file,
    units_override_check,
)
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    requires_ds,
    small_patch_amr,
)


_fields = ("dbz", "thrhopert", "zvort")

## To-Do - Need to know the actual location of test
## data for these tests! Probably need to give a small,
## static NetCDF file near the tornado.
cm1sim = "CM1Tornado/test_dataset.nc"


@requires_ds(cm1sim, big_data=True)
def test_tornado():
    ds = data_dir_load(sloshing)
    assert_equal(str(ds), "test_dataset.nc")
    for test in small_patch_amr(ds, _fields):
        test_tprmadp.__name__ = test.description
        yield test


@requires_file(cm1sim)
def test_CM1Dataset():
    assert isinstance(data_dir_load(cm1sim), CM1Dataset)


@requires_file(cm1sim)
def test_units_override():
    units_override_check(cm1sim)


@requires_file(cm1sim)
def test_tornado_dataset():
    ds = data_dir_load(cm1sim)
    ## To-Do: Static tests for the specific
    ## NetCDF file given above!
    assert_equal(ds.parameters["time"], 751000000000.0)
    assert_equal(ds.domain_dimensions, np.array([8, 8, 8]))
    assert_equal(ds.domain_left_edge, ds.arr([-2e18, -2e18, -2e18], "code_length"))

    assert_equal(ds.index.num_grids, 73)
    dd = ds.all_data()
    dd["density"]


