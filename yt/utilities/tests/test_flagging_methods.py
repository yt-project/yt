from yt.testing import *
from yt.utilities.flagging_methods import flagging_method_registry

def setup():
    global pf
    pf = fake_random_pf(64)
    pf.h

def test_over_density():
    od_flag = flagging_method_registry["overdensity"](0.75) 
    criterion = (pf.h.grids[0]["Density"] > 0.75)
    assert( np.all( od_flag(pf.h.grids[0]) == criterion) )
