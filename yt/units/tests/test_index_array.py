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

from yt.units.index_array import IndexArray

def test_index_array_multiplication():
    vals = np.random.random((100, 3))
    index = IndexArray(vals, input_units=[u.km, u.g, u.s])
    index *= index
    
    assert(index.units[0] == (u.km * u.km).units)
    assert(index.units[1] == (u.g * u.g).units)
    assert(index.units[2] == (u.s * u.s).units)

    index = IndexArray(vals, input_units=[u.km, u.g, u.s])

    index *= u.km

    assert(index.units[0] == (u.km * u.km).units)
    assert(index.units[1] == (u.km * u.g).units)
    assert(index.units[2] == (u.km * u.s).units)
