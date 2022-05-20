from unittest import SkipTest

import matplotlib

from yt.testing import requires_backend

active_backend = matplotlib.get_backend()
inactive_backend = ({"gtkagg", "macosx", "wx", "tkagg"} - {active_backend}).pop()


def test_requires_inactive_backend():
    @requires_backend(inactive_backend)
    def foo():
        return

    try:
        foo()
    except SkipTest:
        pass
    else:
        raise AssertionError(
            "@requires_backend appears to be broken (skip was expected)"
        )


def test_requires_active_backend():
    @requires_backend(active_backend)
    def foo():
        return

    try:
        foo()
    except SkipTest:
        raise AssertionError(
            "@requires_backend appears to be broken (skip was not expected)"
        )
