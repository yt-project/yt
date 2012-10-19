from yt.testing import *
from yt.utilities.lib import obtain_rvec, obtain_rv_vec

_fields = ("Density", "x-velocity", "y-velocity", "z-velocity")

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

    assert_array_equal(vels[0,:], dd['x-velocity'])
    assert_array_equal(vels[1,:], dd['y-velocity'])
    assert_array_equal(vels[2,:], dd['z-velocity'])
