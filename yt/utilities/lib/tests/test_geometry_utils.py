from yt.testing import *
from yt.utilities.lib.api import obtain_rvec, obtain_rv_vec

_fields = ("density", "velocity_x", "velocity_y", "velocity_z")

def test_obtain_rvec():
    pf = fake_random_pf(64, nprocs=8, fields=_fields, 
           negative = [False, True, True, True])
    
    dd = pf.h.sphere((0.5,0.5,0.5), 0.2)

    coords = obtain_rvec(dd)

    r = np.sqrt(np.sum(coords*coords,axis=0))

    assert_array_less(r.max(), 0.2)

    assert_array_less(0.0, r.min())

def test_obtain_rv_vec():
    pf = fake_random_pf(64, nprocs=8, fields=_fields, 
           negative = [False, True, True, True])

    dd = pf.h.all_data()

    vels = obtain_rv_vec(dd)

    assert_array_equal(vels[0,:], dd['velocity_x'])
    assert_array_equal(vels[1,:], dd['velocity_y'])
    assert_array_equal(vels[2,:], dd['velocity_z'])
