from yt.testing import *
from yt.data_objects.profiles import BinnedProfile1D
from numpy.random import uniform

def setup():
    global pf
    pf = fake_random_pf(64, nprocs=8)
    pf.h
    
def test_update_data() :
    dims = (32,32,32)
    grid_data = [{"Temperature":uniform(size=dims)}
                 for i in xrange(pf.h.num_grids)]
    pf.h.update_data(grid_data)
    prj = pf.h.proj(2, "Temperature")
    prj["Temperature"]
    dd = pf.h.all_data()
    profile = BinnedProfile1D(dd, 10, "Density",
                              dd["Density"].min(),
                              dd["Density"].max())
    profile.add_fields(["Temperature"])
    profile["Temperature"]
                              
