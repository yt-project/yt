import numpy as np

from yt.testing import fake_random_ds
from yt.utilities.flagging_methods import flagging_method_registry


def test_over_density():
    ds = fake_random_ds(64)
    ds.index
    od_flag = flagging_method_registry["overdensity"](0.75) 
    criterion = (ds.index.grids[0]["density"] > 0.75)
    assert( np.all( od_flag(ds.index.grids[0]) == criterion) )
