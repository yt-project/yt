import numpy as np

from yt.testing import assert_equal, fake_amr_ds, fake_random_ds


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def test_mean_sum_integrate():
    for nprocs in [-1, 1, 2, 16]:
        if nprocs == -1:
            ds = fake_amr_ds(fields=("density",), units=("g/cm**3",), particles=20)
        else:
            ds = fake_random_ds(
                32, nprocs=nprocs, fields=("density",), units=("g/cm**3",), particles=20
            )
        ad = ds.all_data()

        # Sums
        q = ad.sum(("gas", "density"))

        q1 = ad.quantities.total_quantity(("gas", "density"))

        assert_equal(q, q1)

        q = ad.sum(("all", "particle_ones"))

        q1 = ad.quantities.total_quantity(("all", "particle_ones"))

        assert_equal(q, q1)

        # Weighted Averages
        w = ad.mean(("gas", "density"))

        w1 = ad.quantities.weighted_average_quantity(
            ("gas", "density"), ("index", "ones")
        )

        assert_equal(w, w1)

        w = ad.mean(("gas", "density"), weight=("gas", "density"))

        w1 = ad.quantities.weighted_average_quantity(
            ("gas", "density"), ("gas", "density")
        )

        assert_equal(w, w1)

        w = ad.mean(("all", "particle_mass"))

        w1 = ad.quantities.weighted_average_quantity(
            ("all", "particle_mass"), ("all", "particle_ones")
        )

        assert_equal(w, w1)

        w = ad.mean(("all", "particle_mass"), weight=("all", "particle_mass"))

        w1 = ad.quantities.weighted_average_quantity(
            ("all", "particle_mass"), ("all", "particle_mass")
        )

        assert_equal(w, w1)

        # Projections
        p = ad.sum(("gas", "density"), axis=0)

        p1 = ds.proj(("gas", "density"), 0, data_source=ad, method="sum")

        assert_equal(p[("gas", "density")], p1[("gas", "density")])

        # Check by axis-name
        p = ad.sum(("gas", "density"), axis="x")

        assert_equal(p[("gas", "density")], p1[("gas", "density")])

        # Now we check proper projections
        p = ad.integrate(("gas", "density"), axis=0)
        p1 = ds.proj(("gas", "density"), 0, data_source=ad)

        assert_equal(p[("gas", "density")], p1[("gas", "density")])

        # Check by axis-name
        p = ad.integrate(("gas", "density"), axis="x")

        assert_equal(p[("gas", "density")], p1[("gas", "density")])


def test_min_max():
    for nprocs in [-1, 1, 2, 16]:
        fields = ["density", "temperature"]
        units = ["g/cm**3", "K"]
        if nprocs == -1:
            ds = fake_amr_ds(fields=fields, units=units, particles=20)
        else:
            ds = fake_random_ds(
                32, nprocs=nprocs, fields=fields, units=units, particles=20
            )

        ad = ds.all_data()

        q = ad.min(("gas", "density")).v
        assert_equal(q, ad[("gas", "density")].min())

        q = ad.max(("gas", "density")).v
        assert_equal(q, ad[("gas", "density")].max())

        q = ad.min(("all", "particle_mass")).v
        assert_equal(q, ad[("all", "particle_mass")].min())

        q = ad.max(("all", "particle_mass")).v
        assert_equal(q, ad[("all", "particle_mass")].max())

        ptp = ad.ptp(("gas", "density")).v
        assert_equal(ptp, ad[("gas", "density")].max() - ad[("gas", "density")].min())

        ptp = ad.ptp(("all", "particle_mass")).v
        assert_equal(
            ptp, ad[("all", "particle_mass")].max() - ad[("all", "particle_mass")].min()
        )

        p = ad.max(("gas", "density"), axis=1)
        p1 = ds.proj(("gas", "density"), 1, data_source=ad, method="max")
        assert_equal(p[("gas", "density")], p1[("gas", "density")])

        p = ad.min(("gas", "density"), axis=1)
        p1 = ds.proj(("gas", "density"), 1, data_source=ad, method="min")
        assert_equal(p[("gas", "density")], p1[("gas", "density")])

        p = ad.max(("gas", "density"), axis="y")
        p1 = ds.proj(("gas", "density"), 1, data_source=ad, method="max")
        assert_equal(p[("gas", "density")], p1[("gas", "density")])

        p = ad.min(("gas", "density"), axis="y")
        p1 = ds.proj(("gas", "density"), 1, data_source=ad, method="min")
        assert_equal(p[("gas", "density")], p1[("gas", "density")])

        # Test that we can get multiple in a single pass

        qrho, qtemp = ad.max([("gas", "density"), ("gas", "temperature")])
        assert_equal(qrho, ad[("gas", "density")].max())
        assert_equal(qtemp, ad[("gas", "temperature")].max())

        qrho, qtemp = ad.min([("gas", "density"), ("gas", "temperature")])
        assert_equal(qrho, ad[("gas", "density")].min())
        assert_equal(qtemp, ad[("gas", "temperature")].min())


def test_argmin():
    fields = ["density", "temperature"]
    units = ["g/cm**3", "K"]
    for nprocs in [-1, 1, 2, 16]:
        if nprocs == -1:
            ds = fake_amr_ds(fields=fields, units=units)
        else:
            ds = fake_random_ds(
                32,
                nprocs=nprocs,
                fields=fields,
                units=units,
            )

        ad = ds.all_data()

        q = ad.argmin(("gas", "density"), axis=[("gas", "density")])
        assert_equal(q, ad[("gas", "density")].min())

        q1, q2 = ad.argmin(
            ("gas", "density"), axis=[("gas", "density"), ("gas", "temperature")]
        )
        mi = np.argmin(ad[("gas", "density")])
        assert_equal(q1, ad[("gas", "density")].min())
        assert_equal(q2, ad[("gas", "temperature")][mi])

        pos = ad.argmin(("gas", "density"))
        mi = np.argmin(ad[("gas", "density")])
        assert_equal(pos[0], ad[("index", "x")][mi])
        assert_equal(pos[1], ad[("index", "y")][mi])
        assert_equal(pos[2], ad[("index", "z")][mi])


def test_argmax():
    fields = ["density", "temperature"]
    units = ["g/cm**3", "K"]
    for nprocs in [-1, 1, 2, 16]:
        if nprocs == -1:
            ds = fake_amr_ds(fields=fields, units=units)
        else:
            ds = fake_random_ds(
                32,
                nprocs=nprocs,
                fields=fields,
                units=units,
            )

        ad = ds.all_data()

        q = ad.argmax(("gas", "density"), axis=[("gas", "density")])
        assert_equal(q, ad[("gas", "density")].max())

        q1, q2 = ad.argmax(
            ("gas", "density"), axis=[("gas", "density"), ("gas", "temperature")]
        )
        mi = np.argmax(ad[("gas", "density")])
        assert_equal(q1, ad[("gas", "density")].max())
        assert_equal(q2, ad[("gas", "temperature")][mi])

        pos = ad.argmax(("gas", "density"))
        mi = np.argmax(ad[("gas", "density")])
        assert_equal(pos[0], ad[("index", "x")][mi])
        assert_equal(pos[1], ad[("index", "y")][mi])
        assert_equal(pos[2], ad[("index", "z")][mi])
