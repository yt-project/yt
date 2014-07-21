from yt.testing import *
import numpy as np

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_extrema():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs = nprocs, fields = ("density",
                "velocity_x", "velocity_y", "velocity_z"))
        sp = ds.sphere("c", (0.25, 'unitary'))
        mi, ma = sp.quantities["Extrema"]("density")
        yield assert_equal, mi, np.nanmin(sp["density"])
        yield assert_equal, ma, np.nanmax(sp["density"])
        dd = ds.all_data()
        mi, ma = dd.quantities["Extrema"]("density")
        yield assert_equal, mi, np.nanmin(dd["density"])
        yield assert_equal, ma, np.nanmax(dd["density"])
        sp = ds.sphere("max", (0.25, 'unitary'))
        yield assert_equal, np.any(np.isnan(sp["radial_velocity"])), False
        mi, ma = dd.quantities["Extrema"]("radial_velocity")
        yield assert_equal, mi, np.nanmin(dd["radial_velocity"])
        yield assert_equal, ma, np.nanmax(dd["radial_velocity"])

def test_average():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs = nprocs, fields = ("density",))
        ad = ds.all_data()
        
        my_mean = ad.quantities["WeightedAverageQuantity"]("density", "ones")
        yield assert_rel_equal, my_mean, ad["density"].mean(), 12

        my_mean = ad.quantities["WeightedAverageQuantity"]("density", "cell_mass")
        a_mean = (ad["density"] * ad["cell_mass"]).sum() / ad["cell_mass"].sum()
        yield assert_rel_equal, my_mean, a_mean, 12

def test_variance():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs = nprocs, fields = ("density", ))
        ad = ds.all_data()
        
        my_std, my_mean = ad.quantities["WeightedVariance"]("density", "ones")
        yield assert_rel_equal, my_mean, ad["density"].mean(), 12
        yield assert_rel_equal, my_std, ad["density"].std(), 12

        my_std, my_mean = ad.quantities["WeightedVariance"]("density", "cell_mass")        
        a_mean = (ad["density"] * ad["cell_mass"]).sum() / ad["cell_mass"].sum()
        yield assert_rel_equal, my_mean, a_mean, 12
        a_std = np.sqrt((ad["cell_mass"] * (ad["density"] - a_mean)**2).sum() / 
                        ad["cell_mass"].sum())
        yield assert_rel_equal, my_std, a_std, 12

if __name__ == "__main__":
    for i in test_extrema():
        i[0](*i[1:])
