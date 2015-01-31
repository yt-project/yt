from yt.testing import *
from yt.utilities.flagging_methods import flagging_method_registry

def setup():
    global ds
    ds = fake_random_ds(64)
    ds.index

def test_over_density():
    od_flag = flagging_method_registry["overdensity"](0.75) 
    criterion = (ds.index.grids[0]["density"] > 0.75)
    assert( np.all( od_flag(ds.index.grids[0]) == criterion) )
