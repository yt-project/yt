from yt.testing import *
from yt.data_objects.profiles import \
    BinnedProfile1D, BinnedProfile2D, BinnedProfile3D

_fields = ("Density", "Temperature", "Dinosaurs", "Tribbles")

def test_profiles():
    pf = fake_random_pf(64, nprocs = 8, fields = _fields)
    nv = pf.domain_dimensions.prod()
    dd = pf.h.all_data()
    (rmi, rma), (tmi, tma), (dmi, dma) = dd.quantities["Extrema"](
        ["Density", "Temperature", "Dinosaurs"])
    rt, tt, dt = dd.quantities["TotalQuantity"](
        ["Density", "Temperature", "Dinosaurs"])
    # First we look at the 
    for nb in [8, 16, 32, 64]:
        for lr in [True, False]:
            # We log all the fields or don't log 'em all.  No need to do them
            # individually.
            for lf in [True, False]: 
                # We have the min and the max, but to avoid cutting them off
                # since we aren't doing end-collect, we cut a bit off the edges
                for ec, e1, e2 in [(False, 0.9, 1.1), (True, 1.0, 1.0)]:
                    p1d = BinnedProfile1D(dd, 
                        nb, "Density", rmi*e1, rma*e2, lf,
                        lr, end_collect=ec)
                    p1d.add_fields(["Ones", "Temperature"], weight=None)
                    yield assert_equal, p1d["Ones"].sum(), nv
                    yield assert_rel_equal, tt, p1d["Temperature"].sum(), 7

                    p2d = BinnedProfile2D(dd, 
                        nb, "Density", rmi*e1, rma*e2, lf,
                        nb, "Temperature", tmi*e1, tma*e2, lf,
                        lr, end_collect=ec)
                    p2d.add_fields(["Ones", "Temperature"], weight=None)
                    yield assert_equal, p2d["Ones"].sum(), nv
                    yield assert_rel_equal, tt, p2d["Temperature"].sum(), 7

                    p3d = BinnedProfile3D(dd, 
                        nb, "Density", rmi*e1, rma*e2, lf,
                        nb, "Temperature", tmi*e1, tma*e2, lf,
                        nb, "Dinosaurs", dmi*e1, dma*e2, lf,
                        lr, end_collect=ec)
                    p3d.add_fields(["Ones", "Temperature"], weight=None)
                    yield assert_equal, p3d["Ones"].sum(), nv
                    yield assert_rel_equal, tt, p3d["Temperature"].sum(), 7

            p1d = BinnedProfile1D(dd, nb, "x", 0.0, 1.0, log_space=False)
            p1d.add_fields("Ones", weight=None)
            av = nv / nb
            yield assert_equal, p1d["Ones"][:-1], np.ones(nb)*av
            # We re-bin ones with a weight now
            p1d.add_fields(["Ones"], weight="Temperature")
            yield assert_equal, p1d["Ones"][:-1], np.ones(nb)

            p2d = BinnedProfile2D(dd, nb, "x", 0.0, 1.0, False,
                                      nb, "y", 0.0, 1.0, False)
            p2d.add_fields("Ones", weight=None)
            av = nv / nb**2
            yield assert_equal, p2d["Ones"][:-1,:-1], np.ones((nb, nb))*av
            # We re-bin ones with a weight now
            p2d.add_fields(["Ones"], weight="Temperature")
            yield assert_equal, p2d["Ones"][:-1,:-1], np.ones((nb, nb))

            p3d = BinnedProfile3D(dd, nb, "x", 0.0, 1.0, False,
                                      nb, "y", 0.0, 1.0, False,
                                      nb, "z", 0.0, 1.0, False)
            p3d.add_fields("Ones", weight=None)
            av = nv / nb**3
            yield assert_equal, p3d["Ones"][:-1,:-1,:-1], np.ones((nb, nb, nb))*av
            # We re-bin ones with a weight now
            p3d.add_fields(["Ones"], weight="Temperature")
            yield assert_equal, p3d["Ones"][:-1,:-1,:-1], np.ones((nb,nb,nb))

