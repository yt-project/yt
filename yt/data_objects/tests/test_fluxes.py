from yt.testing import *

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_flux_calculation():
    pf = fake_random_pf(64, nprocs = 4)
    dd = pf.h.all_data()
    surf = pf.h.surface(dd, "x", 0.51)
    yield assert_equal, surf["x"], 0.51
    flux = surf.calculate_flux("Ones", "Zeros", "Zeros", "Ones")
    yield assert_almost_equal, flux, 1.0, 12

