from yt.testing import *
from yt.data_objects.profiles import BinnedProfile1D
from numpy.random import uniform

def test_update_data() :
    ds = fake_random_ds(64, nprocs=8)
    ds.index
    dims = (32,32,32)
    grid_data = [{"temperature":uniform(size=dims)}
                 for i in range(ds.index.num_grids)]
    ds.index.update_data(grid_data, {('gas', 'temperature'):'K'})
    prj = ds.proj("temperature", 2)
    prj["temperature"]
    dd = ds.all_data()
    profile = BinnedProfile1D(dd, 10, "density",
                              dd["density"].min(),
                              dd["density"].max())
    profile.add_fields(["temperature"])
    profile["temperature"]
                              
