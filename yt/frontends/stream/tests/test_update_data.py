import numpy as np
from yt.testing import \
    fake_particle_ds, \
    fake_random_ds
from yt.data_objects.profiles import create_profile

def test_update_data_grid():
    ds = fake_random_ds(64, nprocs=8)
    ds.index
    dims = (32,32,32)
    grid_data = [{"temperature": np.random.uniform(size=dims)}
                 for i in range(ds.index.num_grids)]
    ds.index.update_data(grid_data)
    prj = ds.proj("temperature", 2)
    prj["temperature"]
    dd = ds.all_data()
    profile = create_profile(dd, "density", "temperature", 10)
    profile["temperature"]

def test_update_data_particle():
    npart = 100
    ds = fake_particle_ds(npart=npart)
    part_data = {"temperature": np.random.rand(npart)}
    ds.index.update_data(part_data)
    assert ("io", "temperature") in ds.field_list
    dd = ds.all_data()
    dd["temperature"]
