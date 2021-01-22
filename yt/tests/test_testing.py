import matplotlib
import numpy as np
import pytest

from yt.config import ytcfg
from yt.testing import assert_equal, requires_backend


def test_requires_backend():
    backend = matplotlib.get_backend().lower()
    other_backends = {"gtkagg", "macosx", "wx", "tkagg"} - {backend}

    @requires_backend(other_backends.pop())
    def plot_a():
        return True

    @requires_backend(backend)
    def plot_b():
        return True

    assert_equal(plot_b(), True)
    if not ytcfg.get("yt", "internals", "within_pytest"):
        assert_equal(plot_a(), None)
    else:
        # NOTE: This doesn't actually work. pytest.skip() doesn't actually
        # raise the exception but rather returns control to the function's
        # (test_requires_backend) caller, breaking immediately. As such,
        # this assert_rasies never actually happens. See the comment
        # in the definition of requires_backend for why pytest.skip is used
        np.testing.assert_raises(plot_a(), pytest.skip.Exception)
