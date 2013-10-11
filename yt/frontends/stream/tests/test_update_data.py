from yt.testing import *
from yt.data_objects.profiles import BinnedProfile1D
from numpy.random import uniform

def test_update_data() :
    pf = fake_random_pf(64, nprocs=8)
    pf.h
    dims = (32,32,32)
    grid_data = [{"temperature":uniform(size=dims)}
                 for i in xrange(pf.h.num_grids)]
    pf.h.update_data(grid_data)
    prj = pf.h.proj("temperature", 2)
    prj["temperature"]
    dd = pf.h.all_data()
    profile = BinnedProfile1D(dd, 10, "density",
                              dd["density"].min(),
                              dd["density"].max())
    profile.add_fields(["temperature"])
    profile["temperature"]
                              
