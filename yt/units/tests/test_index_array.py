"""
Test ndarray subclass that handles indexing along dimensions with units.




"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import yt.units as u
import numpy as np

from yt.testing import assert_equal
from yt.units.index_array import IndexArray
from yt.units.yt_array import YTArray

def test_multiplication():
    vals = np.random.random((100, 3))
    index = IndexArray(vals, input_units=[u.km, u.g, u.s])
    index *= index

    assert_equal(index.units[0], (u.km * u.km).units)
    assert_equal(index.units[1], (u.g * u.g).units)
    assert_equal(index.units[2], (u.s * u.s).units)

    index = IndexArray(vals, input_units=[u.km, u.g, u.s])

    index *= u.km

    assert_equal(index.units[0], (u.km * u.km).units)
    assert_equal(index.units[1], (u.km * u.g).units)
    assert_equal(index.units[2], (u.km * u.s).units)

    index2 = index * 2

    assert_equal(index.units[0], index2.units[0])
    assert_equal(index.units[1], index2.units[1])
    assert_equal(index.units[2], index2.units[2])

    index3 = 2 * index

    assert_equal(index.units[0], index3.units[0])
    assert_equal(index.units[1], index3.units[1])
    assert_equal(index.units[2], index3.units[2])


def test_slicing():
    vals = np.arange(300)
    vals.shape = (100, 3)
    index = IndexArray(vals, input_units=[u.km, u.g, u.s])

    ret1 = YTArray(3*np.arange(100), 'km')
    ret2 = IndexArray(np.arange(3), [u.km, u.g, u.s])

    assert_equal(index[:, 0], ret1)
    assert_equal(index[0, :], ret2)
