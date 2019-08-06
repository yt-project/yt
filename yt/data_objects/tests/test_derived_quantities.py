import numpy as np

from yt.testing import \
    fake_random_ds, \
    fake_particle_ds, \
    assert_equal, \
    assert_rel_equal, \
    assert_almost_equal

from yt import particle_filter

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_extrema():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs = nprocs, fields = ("density",
                "velocity_x", "velocity_y", "velocity_z"))
        for sp in [ds.sphere("c", (0.25, 'unitary')), ds.r[0.5,:,:]]:
            mi, ma = sp.quantities["Extrema"]("density")
            assert_equal(mi, np.nanmin(sp["density"]))
            assert_equal(ma, np.nanmax(sp["density"]))
            dd = ds.all_data()
            mi, ma = dd.quantities["Extrema"]("density")
            assert_equal(mi, np.nanmin(dd["density"]))
            assert_equal(ma, np.nanmax(dd["density"]))
            sp = ds.sphere("max", (0.25, 'unitary'))
            assert_equal(np.any(np.isnan(sp["radial_velocity"])), False)
            mi, ma = dd.quantities["Extrema"]("radial_velocity")
            assert_equal(mi, np.nanmin(dd["radial_velocity"]))
            assert_equal(ma, np.nanmax(dd["radial_velocity"]))

def test_average():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs = nprocs, fields = ("density",))
        for ad in [ds.all_data(), ds.r[0.5, :, :]]:

            my_mean = ad.quantities["WeightedAverageQuantity"]("density", "ones")
            assert_rel_equal(my_mean, ad["density"].mean(), 12)

            my_mean = ad.quantities["WeightedAverageQuantity"]("density", "cell_mass")
            a_mean = (ad["density"] * ad["cell_mass"]).sum() / ad["cell_mass"].sum()
            assert_rel_equal(my_mean, a_mean, 12)

def test_variance():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs = nprocs, fields = ("density", ))
        for ad in [ds.all_data(), ds.r[0.5, :, :]]:

            my_std, my_mean = ad.quantities["WeightedVariance"]("density", "ones")
            assert_rel_equal(my_mean, ad["density"].mean(), 12)
            assert_rel_equal(my_std, ad["density"].std(), 12)

            my_std, my_mean = ad.quantities["WeightedVariance"]("density", "cell_mass")
            a_mean = (ad["density"] * ad["cell_mass"]).sum() / ad["cell_mass"].sum()
            assert_rel_equal(my_mean, a_mean, 12)
            a_std = np.sqrt((ad["cell_mass"] * (ad["density"] - a_mean)**2).sum() /
                            ad["cell_mass"].sum())
            assert_rel_equal(my_std, a_std, 12)

def test_max_location():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs = nprocs, fields = ("density", ))
        for ad in [ds.all_data(), ds.r[0.5, :, :]]:

            mv, x, y, z = ad.quantities.max_location(("gas", "density"))

            assert_equal(mv, ad["density"].max())

            mi = np.argmax(ad["density"])

            assert_equal(ad["x"][mi], x)
            assert_equal(ad["y"][mi], y)
            assert_equal(ad["z"][mi], z)

def test_min_location():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs = nprocs, fields = ("density", ))
        for ad in [ds.all_data(), ds.r[0.5, :, :]]:

            mv, x, y, z = ad.quantities.min_location(("gas", "density"))

            assert_equal(mv, ad["density"].min())

            mi = np.argmin(ad["density"])

            assert_equal(ad["x"][mi], x)
            assert_equal(ad["y"][mi], y)
            assert_equal(ad["z"][mi], z)

def test_sample_at_min_field_values():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs = nprocs,
            fields = ("density", "temperature", "velocity_x"))
        for ad in [ds.all_data(), ds.r[0.5, :, :]]:

            mv, temp, vm = ad.quantities.sample_at_min_field_values(
                "density", ["temperature", "velocity_x"])

            assert_equal(mv, ad["density"].min())

            mi = np.argmin(ad["density"])

            assert_equal(ad["temperature"][mi], temp)
            assert_equal(ad["velocity_x"][mi], vm)

def test_sample_at_max_field_values():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs = nprocs,
            fields = ("density", "temperature", "velocity_x"))
        for ad in [ds.all_data(), ds.r[0.5, :, :]]:

            mv, temp, vm = ad.quantities.sample_at_max_field_values(
                "density", ["temperature", "velocity_x"])

            assert_equal(mv, ad["density"].max())

            mi = np.argmax(ad["density"])

            assert_equal(ad["temperature"][mi], temp)
            assert_equal(ad["velocity_x"][mi], vm)

def test_derived_quantities_with_particle_types():

    ds = fake_particle_ds()

    @particle_filter(requires=["particle_position_x"], filtered_type='all')
    def low_x(pfilter,data):
        return data['particle_position_x'].in_units('code_length')<0.5
    ds.add_particle_filter('low_x')

    ad=ds.all_data()

    for ptype in ['all','low_x']:
        #Check bulk velocity
        bulk_vx=(ad[(ptype,'particle_mass')]*ad[(ptype,'particle_velocity_x')]/ad[(ptype,'particle_mass')].sum()).sum()
        assert_almost_equal(ad.quantities.bulk_velocity(use_gas=False,use_particles=True,particle_type=ptype)[0],bulk_vx,5)

        #Check center of mass
        com_x=(ad[(ptype,'particle_mass')]*ad[(ptype,'particle_position_x')]/ad[(ptype,'particle_mass')].sum()).sum()
        assert_almost_equal(ad.quantities.center_of_mass(use_gas=False,use_particles=True,particle_type=ptype)[0],com_x,5)

        #Check angular momentum vector
        l_x=(ad[(ptype,'particle_specific_angular_momentum_x')]*ad[(ptype,'particle_mass')]/ad[(ptype,'particle_mass')].sum()).sum()
        assert_almost_equal(ad.quantities.angular_momentum_vector(use_gas=False,use_particles=True,particle_type=ptype)[0],l_x,5)

    #Check spin parameter values
    assert_almost_equal(ad.quantities.spin_parameter(use_gas=False,use_particles=True),655.7311454765503)
    assert_almost_equal(ad.quantities.spin_parameter(use_gas=False,use_particles=True,particle_type='low_x'),1309.164886405665)
