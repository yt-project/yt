from yt.testing import \
    fake_random_ds, \
    assert_almost_equal, \
    assert_equal

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_flux_calculation():
    ds = fake_random_ds(64, nprocs = 4)
    dd = ds.all_data()
    surf = ds.surface(dd, "x", 0.51)
    yield assert_equal, surf["x"], 0.51
    flux = surf.calculate_flux("ones", "zeros", "zeros", "ones")
    yield assert_almost_equal, flux, 1.0, 12

def test_sampling():
    ds = fake_random_ds(64, nprocs = 4)
    dd = ds.all_data()
    for i, ax in enumerate('xyz'):
        surf = ds.surface(dd, ax, 0.51)
        surf.get_data(ax, "vertex")
        yield assert_equal, surf.vertex_samples[ax], surf.vertices[i,:]
