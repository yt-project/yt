from yt.testing import fake_random_ds, assert_equal
import numpy as np


def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_mean_and_sum():
    for nprocs in [1]:
        ds = fake_random_ds(16, nprocs=nprocs,
                            fields=("density",))
        ad = ds.all_data()

        # Sums
        q = ad.sum('density')

        q1 = ad.quantities.total_quantity('density')

        yield assert_equal, q, q1

        q2 = ad.mean('density', weight=None)

        yield assert_equal, q, q2

        # Weighted Averages
        w = ad.mean("density")

        w1 = ad.quantities.weighted_average_quantity('density', 'ones')

        yield assert_equal, w, w1

        w = ad.mean("density", weight="density")

        w1 = ad.quantities.weighted_average_quantity('density', 'density')

        yield assert_equal, w, w1

        # Projections
        p = ad.sum('density', axis=0)

        p1 = ds.proj('density', 0, data_source=ad)

        yield assert_equal, p['density'], p1['density']


if __name__ == "__main__":
    for args in test_mean_and_sum():
        args[0](*args[1:])
