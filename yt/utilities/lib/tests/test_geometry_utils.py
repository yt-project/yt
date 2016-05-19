import numpy as np

from yt.testing import \
    fake_random_ds, \
    assert_array_less, \
    assert_array_equal
from yt.utilities.lib.misc_utilities import obtain_rvec, obtain_rv_vec

_fields = ("density", "velocity_x", "velocity_y", "velocity_z")

def test_obtain_rvec():
    ds = fake_random_ds(64, nprocs=8, fields=_fields, 
           negative = [False, True, True, True])
    
    dd = ds.sphere((0.5,0.5,0.5), 0.2)

    coords = obtain_rvec(dd)

    r = np.sqrt(np.sum(coords*coords,axis=0))

    assert_array_less(r.max(), 0.2)

    assert_array_less(0.0, r.min())

def test_obtain_rv_vec():
    ds = fake_random_ds(64, nprocs=8, fields=_fields, 
           negative = [False, True, True, True])

    dd = ds.all_data()

    vels = obtain_rv_vec(dd)

    assert_array_equal(vels[0,:], dd['velocity_x'])
    assert_array_equal(vels[1,:], dd['velocity_y'])
    assert_array_equal(vels[2,:], dd['velocity_z'])
