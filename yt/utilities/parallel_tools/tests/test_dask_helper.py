from dask import array, compute

import yt.utilities.parallel_tools.dask_helper as dh
from yt.testing import assert_equal


def test_compute_helper():
    # ensure that results returned by dask_helper.compute match those
    # returned by standard dask compute methods.
    f = array.random.random((1000, 5), chunks=500)
    meanval = f.mean()
    # standard compute methods:
    r_orig = meanval.compute()
    r_orig2 = compute(meanval)[0]

    # wrapped method
    new = dh.compute(meanval)[0]

    assert_equal(r_orig, new)
    assert_equal(r_orig2, new)


def test_client_spin_up_down():
    c = dh.dask_client
    c.start_client()
    assert_equal("running", c.client.status)
    c.close()
    assert_equal(None, c.client)
