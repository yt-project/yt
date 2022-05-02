import matplotlib

from yt.testing import requires_backend

active_backend = matplotlib.get_backend()
inactive_backend = ({"gtkagg", "macosx", "wx", "tkagg"} - {active_backend}).pop()


@requires_backend(inactive_backend)
def test_requires_inactive_backend():
    # this test cannot pass, check that it's always skipped
    raise AssertionError("@requires_backend appears to be broken (expected skip)")


def test_requires_active_backend():
    from unittest import SkipTest

    @requires_backend(active_backend)
    def foo():
        assert True

    try:
        foo()
    except SkipTest:
        raise AssertionError(
            "@requires_backend appears to be broken (skip was not expected)"
        )
