"""
Tests for yt.testing
"""
#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import matplotlib

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

    assert_equal(plot_a(), None)
    assert_equal(plot_b(), True)
