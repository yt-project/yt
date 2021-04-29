import numpy as np

from yt.testing import assert_equal, assert_rel_equal, fake_random_ds


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def test_data_collection():
    # We decompose in different ways
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs=nprocs)
        coll = ds.data_collection(ds.index.grids)
        crho = coll[("gas", "density")].sum(dtype="float64").to_ndarray()
        grho = np.sum(
            [g[("gas", "density")].sum(dtype="float64") for g in ds.index.grids],
            dtype="float64",
        )
        assert_rel_equal(np.array([crho]), np.array([grho]), 12)
        assert_equal(coll.size, ds.domain_dimensions.prod())
        for gi in range(ds.index.num_grids):
            grids = ds.index.grids[: gi + 1]
            coll = ds.data_collection(grids)
            crho = coll[("gas", "density")].sum(dtype="float64")
            grho = np.sum(
                [g[("gas", "density")].sum(dtype="float64") for g in grids],
                dtype="float64",
            )
            assert_rel_equal(np.array([crho]), np.array([grho]), 12)
            assert_equal(coll.size, sum(g.ActiveDimensions.prod() for g in grids))
