import numpy as np

from yt.testing import \
    fake_random_ds, \
    assert_equal, \
    assert_rel_equal

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_data_collection():
    # We decompose in different ways
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(16, nprocs = nprocs)
        coll = ds.data_collection(ds.index.grids)
        crho = coll["density"].sum(dtype="float64").to_ndarray()
        grho = np.sum([g["density"].sum(dtype="float64") for g in ds.index.grids],
                      dtype="float64")
        yield assert_rel_equal, np.array([crho]), np.array([grho]), 12
        yield assert_equal, coll.size, ds.domain_dimensions.prod()
        for gi in range(ds.index.num_grids):
            grids = ds.index.grids[:gi+1]
            coll = ds.data_collection(grids)
            crho = coll["density"].sum(dtype="float64")
            grho = np.sum([g["density"].sum(dtype="float64") for g in grids],
                          dtype="float64")
            yield assert_rel_equal, np.array([crho]), np.array([grho]), 12
            yield assert_equal, coll.size, \
                    sum(g.ActiveDimensions.prod() for g in grids)
