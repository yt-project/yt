from yt.testing import fake_random_ds, fake_amr_ds, assert_equal
import numpy as np


def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_mean_sum_integrate():
    for nprocs in [-1, 1, 2, 16]:
        if nprocs == -1:
            ds = fake_amr_ds(fields=("density",), particles=20)
        else:
            ds = fake_random_ds(32, nprocs=nprocs, fields=("density",),
                                particles=20)
        ad = ds.all_data()

        # Sums
        q = ad.sum('density')

        q1 = ad.quantities.total_quantity('density')

        assert_equal(q, q1)

        q = ad.sum('particle_ones')

        q1 = ad.quantities.total_quantity('particle_ones')

        assert_equal(q, q1)

        # Weighted Averages
        w = ad.mean("density")

        w1 = ad.quantities.weighted_average_quantity('density', 'ones')

        assert_equal(w, w1)

        w = ad.mean("density", weight="density")

        w1 = ad.quantities.weighted_average_quantity('density', 'density')

        assert_equal(w, w1)

        w = ad.mean('particle_mass')

        w1 = ad.quantities.weighted_average_quantity(
            'particle_mass', 'particle_ones')

        assert_equal(w, w1)

        w = ad.mean('particle_mass', weight='particle_mass')

        w1 = ad.quantities.weighted_average_quantity(
            'particle_mass', 'particle_mass')

        assert_equal(w, w1)

        # Projections
        p = ad.sum('density', axis=0)

        p1 = ds.proj('density', 0, data_source=ad, method="sum")

        assert_equal(p['density'], p1['density'])

        # Check by axis-name
        p = ad.sum('density', axis='x')

        assert_equal(p['density'], p1['density'])

        # Now we check proper projections
        p = ad.integrate("density", axis=0)
        p1 = ds.proj("density", 0, data_source=ad)

        assert_equal(p['density'], p1['density'])

        # Check by axis-name
        p = ad.integrate('density', axis='x')

        assert_equal(p['density'], p1['density'])

def test_min_max():
    for nprocs in [-1, 1, 2, 16]:
        if nprocs == -1:
            ds = fake_amr_ds(fields=("density","temperature"), particles=20)
        else:
            ds = fake_random_ds(32, nprocs=nprocs,
                fields=("density","temperature"), particles=20)

        ad = ds.all_data()

        q = ad.min("density").v
        assert_equal(q, ad["density"].min())

        q = ad.max("density").v
        assert_equal(q, ad["density"].max())

        q = ad.min('particle_mass').v
        assert_equal(q, ad['particle_mass'].min())

        q = ad.max('particle_mass').v
        assert_equal(q, ad['particle_mass'].max())

        ptp = ad.ptp("density").v
        assert_equal(ptp, ad["density"].max() - ad["density"].min())

        ptp = ad.ptp("particle_mass").v
        assert_equal(ptp, ad["particle_mass"].max() - ad["particle_mass"].min())

        p = ad.max("density", axis=1)
        p1 = ds.proj("density", 1, data_source=ad, method="mip")
        assert_equal(p["density"], p1["density"])

        p = ad.max("density", axis="y")
        p1 = ds.proj("density", 1, data_source=ad, method="mip")
        assert_equal(p["density"], p1["density"])

        # Test that we can get multiple in a single pass

        qrho, qtemp = ad.max(["density", "temperature"])
        assert_equal(qrho, ad["density"].max())
        assert_equal(qtemp, ad["temperature"].max())

        qrho, qtemp = ad.min(["density", "temperature"])
        assert_equal(qrho, ad["density"].min())
        assert_equal(qtemp, ad["temperature"].min())

def test_argmin():
    for nprocs in [-1, 1, 2, 16]:
        if nprocs == -1:
            ds = fake_amr_ds(fields=("density","temperature"))
        else:
            ds = fake_random_ds(32, nprocs=nprocs,
                fields=("density","temperature"))

        ad = ds.all_data()

        q = ad.argmin("density", axis=["density"])
        assert_equal(q, ad["density"].min())

        q1, q2 = ad.argmin("density", axis=["density", "temperature"])
        mi = np.argmin(ad["density"])
        assert_equal(q1, ad["density"].min())
        assert_equal(q2, ad["temperature"][mi])

        pos = ad.argmin("density")
        mi = np.argmin(ad["density"])
        assert_equal(pos[0], ad["x"][mi])
        assert_equal(pos[1], ad["y"][mi])
        assert_equal(pos[2], ad["z"][mi])

def test_argmax():
    for nprocs in [-1, 1, 2, 16]:
        if nprocs == -1:
            ds = fake_amr_ds(fields=("density","temperature"))
        else:
            ds = fake_random_ds(32, nprocs=nprocs,
                fields=("density","temperature"))

        ad = ds.all_data()

        q = ad.argmax("density", axis=["density"])
        assert_equal(q, ad["density"].max())

        q1, q2 = ad.argmax("density", axis=["density", "temperature"])
        mi = np.argmax(ad["density"])
        assert_equal(q1, ad["density"].max())
        assert_equal(q2, ad["temperature"][mi])

        pos = ad.argmax("density")
        mi = np.argmax(ad["density"])
        assert_equal(pos[0], ad["x"][mi])
        assert_equal(pos[1], ad["y"][mi])
        assert_equal(pos[2], ad["z"][mi])
