#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.testing import \
    assert_allclose_units, \
    fake_random_ds, \
    requires_file

def test_AM_value():
    ds = fake_random_ds(16)

    sp = ds.sphere([.5]*3, (0.1, 'code_length'))

    x0 = sp.center
    v0 = ds.arr([1, 2, 3], 'km/s')

    sp.set_field_parameter('bulk_velocity', v0)

    X = (ds.arr([sp[k] for k in 'xyz']) - x0[:, None]).T
    V = (ds.arr([sp['velocity_'+k] for k in 'xyz']) - v0[:, None]).T

    sAM_manual = ds.arr(np.cross(V, X), X.units*V.units)
    sAM = ds.arr([sp['specific_angular_momentum_'+k] for k in 'xyz']).T

    assert_allclose_units(sAM_manual, sAM)


