from yt.testing import *

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_data_collection():
    # We decompose in different ways
    for nprocs in [1, 2, 4, 8]:
        pf = fake_random_pf(16, nprocs = nprocs)
        coll = pf.h.data_collection(pf.domain_center, pf.h.grids)
        crho = coll["Density"].sum(dtype="float64")
        grho = np.sum([g["Density"].sum(dtype="float64") for g in pf.h.grids],
                      dtype="float64")
        yield assert_rel_equal, crho, grho, 12
        yield assert_equal, coll.size, pf.domain_dimensions.prod()
        for gi in range(pf.h.num_grids):
            grids = pf.h.grids[:gi+1]
            coll = pf.h.data_collection(pf.domain_center, grids)
            crho = coll["Density"].sum(dtype="float64")
            grho = np.sum([g["Density"].sum(dtype="float64") for g in grids],
                          dtype="float64")
            yield assert_rel_equal, crho, grho, 12
            yield assert_equal, coll.size, \
                    sum(g.ActiveDimensions.prod() for g in grids)
