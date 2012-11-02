from yt.testing import *
from yt.data_objects.profiles import \
    BinnedProfile1D, BinnedProfile2D, BinnedProfile3D

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_covering_grid():
    # We decompose in different ways
    for level in [0, 1, 2]:
        for nprocs in [1, 2, 4, 8]:
            pf = fake_random_pf(16, nprocs = nprocs)
            dn = pf.refine_by**level 
            cg = pf.h.covering_grid(level, [0.0, 0.0, 0.0],
                    dn * pf.domain_dimensions)
            yield assert_equal, cg["Ones"].max(), 1.0
            yield assert_equal, cg["Ones"].min(), 1.0
            yield assert_equal, cg["CellVolume"].sum(), pf.domain_width.prod()
            for g in pf.h.grids:
                di = g.get_global_startindex()
                dd = g.ActiveDimensions
                for i in range(dn):
                    f = cg["Density"][dn*di[0]+i:dn*(di[0]+dd[0])+i:dn,
                                      dn*di[1]+i:dn*(di[1]+dd[1])+i:dn,
                                      dn*di[2]+i:dn*(di[2]+dd[2])+i:dn]
                    yield assert_equal, f, g["Density"]
