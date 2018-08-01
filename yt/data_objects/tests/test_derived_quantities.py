import numpy as np

from yt.testing import \
    fake_random_ds, \
    assert_equal, \
    assert_rel_equal, \
    assert_almost_equal

from yt import particle_filter
from yt import load_particles
from yt.units import parsec, Msun

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

    #ds=fake_particle_data()
    #ad=ds.all_data()

    @particle_filter(requires=["particle_position_x"], filtered_type='all')
    def low_x(pfilter,data):
        return data['particle_position_x'].in_units('code_length')<0.5

    #Create pre-defined particle dataset
    nparts=100
    ppx = pvx = np.array(np.linspace(0,1,nparts))
    ppy = pvy = np.array(np.linspace(1,0,nparts))
    ppz = pvz = ppx-ppy
    pm = np.ones(nparts)

    data={'particle_position_x':ppx,
          'particle_position_y':ppy,
          'particle_position_z':ppx,
          'particle_velocity_x':pvx,
          'particle_velocity_y':pvy,
          'particle_velocity_z':pvz,
          'particle_mass':pm}

    bbox=np.array([[min(ppx), max(ppx)], [min(ppy), max(ppy)], [min(ppz), max(ppz)]])
    ds = load_particles(data, n_ref=6, bbox=bbox, length_unit=parsec, mass_unit=Msun)

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


    #Check spin parameter values
    assert_almost_equal(ad.quantities.spin_parameter(use_gas=False,use_particles=True),1.1973910384754621e+27)
    assert_almost_equal(ad.quantities.spin_parameter(use_gas=False,use_particles=True,particle_type='low_x'),2.1171960652141979e+27)
