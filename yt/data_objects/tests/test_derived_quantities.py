import numpy as np

import yt
from yt import particle_filter
from yt.testing import (
    assert_almost_equal,
    assert_equal,
    assert_rel_equal,
    fake_particle_ds,
    fake_random_ds,
    fake_sph_orientation_ds,
    requires_file,
)


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def test_extrema():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(
            16,
            nprocs=nprocs,
            fields=("density", "velocity_x", "velocity_y", "velocity_z"),
            units=("g/cm**3", "cm/s", "cm/s", "cm/s"),
        )
        for sp in [ds.sphere("c", (0.25, "unitary")), ds.r[0.5, :, :]]:
            mi, ma = sp.quantities["Extrema"](("gas", "density"))
            assert_equal(mi, np.nanmin(sp[("gas", "density")]))
            assert_equal(ma, np.nanmax(sp[("gas", "density")]))
            dd = ds.all_data()
            mi, ma = dd.quantities["Extrema"](("gas", "density"))
            assert_equal(mi, np.nanmin(dd[("gas", "density")]))
            assert_equal(ma, np.nanmax(dd[("gas", "density")]))
            sp = ds.sphere("max", (0.25, "unitary"))
            assert_equal(np.any(np.isnan(sp[("gas", "radial_velocity")])), False)
            mi, ma = dd.quantities["Extrema"](("gas", "radial_velocity"))
            assert_equal(mi, np.nanmin(dd[("gas", "radial_velocity")]))
            assert_equal(ma, np.nanmax(dd[("gas", "radial_velocity")]))


def test_average():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs=nprocs, fields=("density",), units=("g/cm**3",))
        for ad in [ds.all_data(), ds.r[0.5, :, :]]:

            my_mean = ad.quantities["WeightedAverageQuantity"](
                ("gas", "density"), ("index", "ones")
            )
            assert_rel_equal(my_mean, ad[("gas", "density")].mean(), 12)

            my_mean = ad.quantities["WeightedAverageQuantity"](
                ("gas", "density"), ("gas", "cell_mass")
            )
            a_mean = (ad[("gas", "density")] * ad[("gas", "cell_mass")]).sum() / ad[
                ("gas", "cell_mass")
            ].sum()
            assert_rel_equal(my_mean, a_mean, 12)


def test_standard_deviation():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs=nprocs, fields=("density",), units=("g/cm**3",))
        for ad in [ds.all_data(), ds.r[0.5, :, :]]:

            my_std, my_mean = ad.quantities["WeightedStandardDeviation"](
                ("gas", "density"), ("index", "ones")
            )
            assert_rel_equal(my_mean, ad[("gas", "density")].mean(), 12)
            assert_rel_equal(my_std, ad[("gas", "density")].std(), 12)

            my_std, my_mean = ad.quantities["WeightedStandardDeviation"](
                ("gas", "density"), ("gas", "cell_mass")
            )
            a_mean = (ad[("gas", "density")] * ad[("gas", "cell_mass")]).sum() / ad[
                ("gas", "cell_mass")
            ].sum()
            assert_rel_equal(my_mean, a_mean, 12)
            a_std = np.sqrt(
                (
                    ad[("gas", "cell_mass")] * (ad[("gas", "density")] - a_mean) ** 2
                ).sum()
                / ad[("gas", "cell_mass")].sum()
            )
            assert_rel_equal(my_std, a_std, 12)


def test_max_location():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs=nprocs, fields=("density",), units=("g/cm**3",))
        for ad in [ds.all_data(), ds.r[0.5, :, :]]:

            mv, x, y, z = ad.quantities.max_location(("gas", "density"))

            assert_equal(mv, ad[("gas", "density")].max())

            mi = np.argmax(ad[("gas", "density")])

            assert_equal(ad[("index", "x")][mi], x)
            assert_equal(ad[("index", "y")][mi], y)
            assert_equal(ad[("index", "z")][mi], z)


def test_min_location():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs=nprocs, fields=("density",), units=("g/cm**3",))
        for ad in [ds.all_data(), ds.r[0.5, :, :]]:

            mv, x, y, z = ad.quantities.min_location(("gas", "density"))

            assert_equal(mv, ad[("gas", "density")].min())

            mi = np.argmin(ad[("gas", "density")])

            assert_equal(ad[("index", "x")][mi], x)
            assert_equal(ad[("index", "y")][mi], y)
            assert_equal(ad[("index", "z")][mi], z)


def test_sample_at_min_field_values():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(
            16,
            nprocs=nprocs,
            fields=("density", "temperature", "velocity_x"),
            units=("g/cm**3", "K", "cm/s"),
        )
        for ad in [ds.all_data(), ds.r[0.5, :, :]]:

            mv, temp, vm = ad.quantities.sample_at_min_field_values(
                ("gas", "density"), [("gas", "temperature"), ("gas", "velocity_x")]
            )

            assert_equal(mv, ad[("gas", "density")].min())

            mi = np.argmin(ad[("gas", "density")])

            assert_equal(ad[("gas", "temperature")][mi], temp)
            assert_equal(ad[("gas", "velocity_x")][mi], vm)


def test_sample_at_max_field_values():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(
            16,
            nprocs=nprocs,
            fields=("density", "temperature", "velocity_x"),
            units=("g/cm**3", "K", "cm/s"),
        )
        for ad in [ds.all_data(), ds.r[0.5, :, :]]:

            mv, temp, vm = ad.quantities.sample_at_max_field_values(
                ("gas", "density"), [("gas", "temperature"), ("gas", "velocity_x")]
            )

            assert_equal(mv, ad[("gas", "density")].max())

            mi = np.argmax(ad[("gas", "density")])

            assert_equal(ad[("gas", "temperature")][mi], temp)
            assert_equal(ad[("gas", "velocity_x")][mi], vm)


def test_in_memory_sph_derived_quantities():
    ds = fake_sph_orientation_ds()
    ad = ds.all_data()

    ang_mom = ad.quantities.angular_momentum_vector()
    assert_equal(ang_mom, [0, 0, 0])

    bv = ad.quantities.bulk_velocity()
    assert_equal(bv, [0, 0, 0])

    com = ad.quantities.center_of_mass()
    assert_equal(com, [1 / 7, (1 + 2) / 7, (1 + 2 + 3) / 7])

    ex = ad.quantities.extrema([("io", "x"), ("io", "y"), ("io", "z")])
    for fex, ans in zip(ex, [[0, 1], [0, 2], [0, 3]]):
        assert_equal(fex, ans)

    for d, v, l in zip("xyz", [1, 2, 3], [[1, 0, 0], [0, 2, 0], [0, 0, 3]]):
        max_d, x, y, z = ad.quantities.max_location(("io", d))
        assert_equal(max_d, v)
        assert_equal([x, y, z], l)

    for d in "xyz":
        min_d, x, y, z = ad.quantities.min_location(("io", d))
        assert_equal(min_d, 0)
        assert_equal([x, y, z], [0, 0, 0])

    tot_m = ad.quantities.total_mass()
    assert_equal(tot_m, [7, 0])

    weighted_av_z = ad.quantities.weighted_average_quantity(("io", "z"), ("io", "z"))
    assert_equal(weighted_av_z, 7 / 3)


iso_collapse = "IsothermalCollapse/snap_505"
tipsy_gal = "TipsyGalaxy/galaxy.00300"


@requires_file(iso_collapse)
@requires_file(tipsy_gal)
def test_sph_datasets_derived_quantities():
    for fname in [tipsy_gal, iso_collapse]:
        ds = yt.load(fname)
        ad = ds.all_data()
        use_particles = "nbody" in ds.particle_types
        ad.quantities.angular_momentum_vector()
        ad.quantities.bulk_velocity(True, use_particles)
        ad.quantities.center_of_mass(True, use_particles)
        ad.quantities.extrema([("gas", "density"), ("gas", "temperature")])
        ad.quantities.min_location(("gas", "density"))
        ad.quantities.max_location(("gas", "density"))
        ad.quantities.total_mass()
        ad.quantities.weighted_average_quantity(("gas", "density"), ("gas", "mass"))


def test_derived_quantities_with_particle_types():

    ds = fake_particle_ds()

    @particle_filter(requires=["particle_position_x"], filtered_type="all")
    def low_x(pfilter, data):
        return (
            data[(pfilter.filtered_type, "particle_position_x")].in_units("code_length")
            < 0.5
        )

    ds.add_particle_filter("low_x")

    ad = ds.all_data()

    for ptype in ["all", "low_x"]:
        # Check bulk velocity
        bulk_vx = (
            ad[(ptype, "particle_mass")]
            * ad[(ptype, "particle_velocity_x")]
            / ad[(ptype, "particle_mass")].sum()
        ).sum()
        assert_almost_equal(
            ad.quantities.bulk_velocity(
                use_gas=False, use_particles=True, particle_type=ptype
            )[0],
            bulk_vx,
            5,
        )

        # Check center of mass
        com_x = (
            ad[(ptype, "particle_mass")]
            * ad[(ptype, "particle_position_x")]
            / ad[(ptype, "particle_mass")].sum()
        ).sum()
        assert_almost_equal(
            ad.quantities.center_of_mass(
                use_gas=False, use_particles=True, particle_type=ptype
            )[0],
            com_x,
            5,
        )

        # Check angular momentum vector
        l_x = (
            ad[(ptype, "particle_specific_angular_momentum_x")]
            * ad[(ptype, "particle_mass")]
            / ad[(ptype, "particle_mass")].sum()
        ).sum()
        assert_almost_equal(
            ad.quantities.angular_momentum_vector(
                use_gas=False, use_particles=True, particle_type=ptype
            )[0],
            l_x,
            5,
        )

    # Check spin parameter values
    assert_almost_equal(
        ad.quantities.spin_parameter(use_gas=False, use_particles=True),
        655.7311454765503,
    )
    assert_almost_equal(
        ad.quantities.spin_parameter(
            use_gas=False, use_particles=True, particle_type="low_x"
        ),
        1309.164886405665,
    )
