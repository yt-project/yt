from yt.testing import *
import numpy as np

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_extrema():
    for nprocs in [1, 2, 4, 8]:
        pf = fake_random_pf(16, nprocs = nprocs, fields = ("Density",
                "x-velocity", "y-velocity", "z-velocity"))
        sp = pf.h.sphere("c", (0.25, '1'))
        (mi, ma), = sp.quantities["Extrema"]("Density")
        yield assert_equal, mi, np.nanmin(sp["Density"])
        yield assert_equal, ma, np.nanmax(sp["Density"])
        dd = pf.h.all_data()
        (mi, ma), = dd.quantities["Extrema"]("Density")
        yield assert_equal, mi, np.nanmin(dd["Density"])
        yield assert_equal, ma, np.nanmax(dd["Density"])
        sp = pf.h.sphere("max", (0.25, '1'))
        yield assert_equal, np.any(np.isnan(sp["RadialVelocity"])), True
        (mi, ma), = dd.quantities["Extrema"]("RadialVelocity")
        yield assert_equal, mi, np.nanmin(dd["RadialVelocity"])
        yield assert_equal, ma, np.nanmax(dd["RadialVelocity"])

def test_average():
    for nprocs in [1, 2, 4, 8]:
        pf = fake_random_pf(16, nprocs = nprocs, fields = ("Density",))
        ad = pf.h.all_data()
        
        my_mean = ad.quantities["WeightedAverageQuantity"]("Density", "ones")
        yield assert_rel_equal, my_mean, ad["Density"].mean(), 12

        my_mean = ad.quantities["WeightedAverageQuantity"]("Density", "CellMass")
        a_mean = (ad["Density"] * ad["CellMass"]).sum() / ad["CellMass"].sum()
        yield assert_rel_equal, my_mean, a_mean, 12

def test_variance():
    for nprocs in [1, 2, 4, 8]:
        pf = fake_random_pf(16, nprocs = nprocs, fields = ("Density", ))
        ad = pf.h.all_data()
        
        my_std, my_mean = ad.quantities["WeightedVariance"]("Density", "ones")
        yield assert_rel_equal, my_mean, ad["Density"].mean(), 12
        yield assert_rel_equal, my_std, ad["Density"].std(), 12

        my_std, my_mean = ad.quantities["WeightedVariance"]("Density", "CellMass")        
        a_mean = (ad["Density"] * ad["CellMass"]).sum() / ad["CellMass"].sum()
        yield assert_rel_equal, my_mean, a_mean, 12
        a_std = np.sqrt((ad["CellMass"] * (ad["Density"] - a_mean)**2).sum() / 
                        ad["CellMass"].sum())
        yield assert_rel_equal, my_std, a_std, 12
