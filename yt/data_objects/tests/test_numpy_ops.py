from yt.testing import fake_random_ds, fake_amr_ds, assert_equal
import numpy as np


def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_mean_sum_integrate():
    for nprocs in [-1, 1, 2, 16]:
        if nprocs == -1:
            ds = fake_amr_ds(fields=("density",))
        else:
            ds = fake_random_ds(32, nprocs=nprocs, fields=("density",))
        ad = ds.all_data()

        # Sums
        q = ad.sum('density')

        q1 = ad.quantities.total_quantity('density')

        yield assert_equal, q, q1

        # Weighted Averages
        w = ad.mean("density")

        w1 = ad.quantities.weighted_average_quantity('density', 'ones')

        yield assert_equal, w, w1

        w = ad.mean("density", weight="density")

        w1 = ad.quantities.weighted_average_quantity('density', 'density')

        yield assert_equal, w, w1

        # Projections
        p = ad.sum('density', axis=0)

        p1 = ds.proj('density', 0, data_source=ad, method="sum")

        yield assert_equal, p['density'], p1['density']

        # Check by axis-name
        p = ad.sum('density', axis='x')

        yield assert_equal, p['density'], p1['density']

        # Now we check proper projections
        p = ad.integrate("density", axis=0)
        p1 = ds.proj("density", 0, data_source=ad)

        yield assert_equal, p['density'], p1['density']

        # Check by axis-name
        p = ad.integrate('density', axis='x')

        yield assert_equal, p['density'], p1['density']

def test_min_max():
    for nprocs in [-1, 1, 2, 16]:
        if nprocs == -1:
            ds = fake_amr_ds(fields=("density","temperature"))
        else:
            ds = fake_random_ds(32, nprocs=nprocs,
                fields=("density","temperature"))

        ad = ds.all_data()

        q = ad.min("density").v
        yield assert_equal, q, ad["density"].min()

        q = ad.max("density").v
        yield assert_equal, q, ad["density"].max()

        ptp = ad.ptp("density").v
        yield assert_equal, ptp, ad["density"].max() - ad["density"].min()

        p = ad.max("density", axis=1)
        p1 = ds.proj("density", 1, data_source=ad, method="mip")
        yield assert_equal, p["density"], p1["density"]

        p = ad.max("density", axis="y")
        p1 = ds.proj("density", 1, data_source=ad, method="mip")
        yield assert_equal, p["density"], p1["density"]

        # Test that we can get multiple in a single pass

        qrho, qtemp = ad.max(["density", "temperature"])
        yield assert_equal, qrho, ad["density"].max()
        yield assert_equal, qtemp, ad["temperature"].max()

        qrho, qtemp = ad.min(["density", "temperature"])
        yield assert_equal, qrho, ad["density"].min()
        yield assert_equal, qtemp, ad["temperature"].min()

def test_argmin():
    for nprocs in [-1, 1, 2, 16]:
        if nprocs == -1:
            ds = fake_amr_ds(fields=("density","temperature"))
        else:
            ds = fake_random_ds(32, nprocs=nprocs,
                fields=("density","temperature"))

        ad = ds.all_data()

        q = ad.argmin("density", axis=["density"])
        yield assert_equal, q, ad["density"].min()

        q1, q2 = ad.argmin("density", axis=["density", "temperature"])
        mi = np.argmin(ad["density"])
        yield assert_equal, q1, ad["density"].min()
        yield assert_equal, q2, ad["temperature"][mi]

        pos = ad.argmin("density")
        mi = np.argmin(ad["density"])
        yield assert_equal, pos[0], ad["x"][mi]
        yield assert_equal, pos[1], ad["y"][mi]
        yield assert_equal, pos[2], ad["z"][mi]

def test_argmax():
    for nprocs in [-1, 1, 2, 16]:
        if nprocs == -1:
            ds = fake_amr_ds(fields=("density","temperature"))
        else:
            ds = fake_random_ds(32, nprocs=nprocs,
                fields=("density","temperature"))

        ad = ds.all_data()

        q = ad.argmax("density", axis=["density"])
        yield assert_equal, q, ad["density"].max()

        q1, q2 = ad.argmax("density", axis=["density", "temperature"])
        mi = np.argmax(ad["density"])
        yield assert_equal, q1, ad["density"].max()
        yield assert_equal, q2, ad["temperature"][mi]

        pos = ad.argmax("density")
        mi = np.argmax(ad["density"])
        yield assert_equal, pos[0], ad["x"][mi]
        yield assert_equal, pos[1], ad["y"][mi]
        yield assert_equal, pos[2], ad["z"][mi]
