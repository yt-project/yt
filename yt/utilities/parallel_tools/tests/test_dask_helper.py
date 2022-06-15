import numpy as np
import pytest

from yt.config import ytcfg
from yt.testing import assert_array_equal, assert_equal
from yt.utilities.on_demand_imports import _dask as dask
from yt.utilities.parallel_tools import dask_helper as dh


@pytest.fixture(autouse=True)
def run_before_and_after_tests():
    # Fixture to execute set dask_enabled config to standard value then reset.
    # ensures that these tests will work even if default value for dask_enabled
    # changes while ensuring other tests will retain the default value even
    # if a test here changes the config value.

    # Setup:
    original_config = ytcfg.get("yt", "internals", "dask_enabled")
    ytcfg.set("yt", "internals", "dask_enabled", False)

    yield  # test will run here

    # Teardown:
    ytcfg.set("yt", "internals", "dask_enabled", original_config)


def test_compute_helper():
    ytcfg.set("yt", "internals", "dask_enabled", True)  # make sure its enabled
    # ensure that results returned by dask_helper.compute match those
    # returned by standard dask compute methods.
    f = dask.array.random.random((1000, 5), chunks=500)
    meanval = f.mean()
    # standard compute methods:
    r_orig = meanval.compute()
    r_orig2 = dask.compute(meanval)[0]

    # wrapped method
    new = dh.compute(meanval)[0]

    assert_equal(r_orig, new)
    assert_equal(r_orig2, new)


def test_client_spin_up_down():
    c = dh.dask_client
    # first check that we catch dask_enabled=False
    with pytest.raises(RuntimeError):
        c.start_client()

    # check that updating dask_enabled works
    ytcfg.set("yt", "internals", "dask_enabled", True)
    assert c._dask_enabled

    # now starting should work
    c.start_client()
    assert_equal("running", c.client.status)
    c.close()
    assert_equal(None, c.client)


def test_dask_config_options():

    # when not enabled, helper functions will pass through arguments or calls
    def delay_func(vals_list):
        results = [vals.min() for vals in vals_list]
        return np.array(results)

    x = [np.linspace(1, 10, 5), np.linspace(2, 30, 5)]
    expected = np.array((1.0, 2.0))
    actual = dh.dask_delay_wrapper(delay_func)(x)  # will compute immediately
    assert_array_equal(expected, actual)

    result = dh.compute(*x)  # will return arguments immediately
    for expected, actual in zip(x, result):
        assert_array_equal(expected, actual)

    # the dh.compute should match dask.compute for non-delayed values
    normal_compute = dask.compute(*x)
    for expected, actual in zip(normal_compute, result):
        assert_array_equal(expected, actual)
