import numpy as np
import pytest

from yt.config import ytcfg
from yt.testing import assert_array_equal, assert_equal
from yt.utilities.on_demand_imports import _dask as dask
from yt.utilities.parallel_tools import dask_helper


@pytest.fixture(autouse=True, scope="module")
def with_dask_enabled():
    # Fixture to set dask_enabled config to standard value then reset.
    # ensures that these tests will work even if default value for dask_enabled
    # changes while ensuring the config value gets reset to the default value
    # after a test runs.

    original_config = ytcfg.get("yt", "internals", "dask_enabled")
    ytcfg.set("yt", "internals", "dask_enabled", False)

    yield  # test will run here

    ytcfg.set("yt", "internals", "dask_enabled", original_config)


@pytest.mark.skipif(not dask.__is_available__, reason="requires dask")
def test_compute_helper(dask_client_fixture):
    ytcfg.set("yt", "internals", "dask_enabled", True)  # make sure its enabled
    # ensure that results returned by dask_helper.compute match those
    # returned by standard dask compute methods.
    f = dask.array.random.random((1000, 5), chunks=500)
    meanval = f.mean()
    # standard compute methods:
    r_orig = meanval.compute()
    r_orig2 = dask.compute(meanval)[0]

    # check the wrapped version
    dask_helper_compute = dask_helper._get_dask_compute(True)
    new = dask_helper_compute(meanval)[0]
    assert_equal(r_orig, new)
    assert_equal(r_orig2, new)

    # when passing compute a non-delayed object, compute just returns the *args,
    # ignoring any kwargs. This checks that we get back compute from dask_helper
    # when expected and also tests that the passthrough function returns the
    # *args just like compute.
    expected = dask.compute(1, 2, 3, this_will_be_dropped=2)
    dask_helper_compute = dask_helper._get_dask_compute(True)
    assert dask_helper_compute == dask.compute
    actual = dask_helper_compute(1, 2, 3, this_will_be_dropped=2)
    assert actual == expected
    dask_helper_compute = dask_helper._get_dask_compute(False)
    assert dask_helper_compute == dask_helper._passthrough_function
    actual = dask_helper_compute(1, 2, 3, this_will_be_dropped=2)
    assert actual == expected

    # check that if dask is not enabled we get back the passthrough function
    ytcfg.set("yt", "internals", "dask_enabled", False)
    dask_helper_compute = dask_helper._get_dask_compute(True)
    assert dask_helper_compute == dask_helper._passthrough_function


@pytest.mark.skipif(not dask.__is_available__, reason="requires dask")
def test_delay_helper(dask_client_fixture):

    # a function to delay (or not)
    def delay_func(vals_list):
        results = [vals.min() for vals in vals_list]
        return np.array(results)

    x = [np.linspace(1, 10, 5), np.linspace(2, 30, 5)]
    expected = np.array((1.0, 2.0))

    # when dask is not enabled, get an empty decorator that calls the function
    # immediately
    delay = dask_helper._get_dask_delayed(True)
    assert delay == dask_helper._empty_decorator
    actual = delay(delay_func)(x)
    assert_array_equal(expected, actual)

    # check that the right functions are selected
    ytcfg.set("yt", "internals", "dask_enabled", True)
    delay = dask_helper._get_dask_delayed(False)
    assert delay == dask_helper._empty_decorator
    delay = dask_helper._get_dask_delayed(True)
    assert delay == dask.delayed

    # check that a delayed operation behaves as expected
    x = [
        dask.array.random.random((1000,), chunks=1000),
        dask.array.random.random((1000,), chunks=1000),
    ]
    assert dask_helper.is_dask_array(x[0])
    expected = dask.delayed(delay_func)(x).compute()
    actual = delay(delay_func)(x).compute()
    assert_array_equal(expected, actual)
