from yt.testing import *
from yt.data_objects.profiles import \
    BinnedProfile1D, BinnedProfile2D, BinnedProfile3D

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_projection():
    for nprocs in [8, 1]:
        # We want to test both 1 proc and 8 procs, to make sure that
        # parallelism isn't broken
        pf = fake_random_pf(64, nprocs = nprocs)
        dims = pf.domain_dimensions
        xn, yn, zn = pf.domain_dimensions
        xi, yi, zi = pf.domain_left_edge + 1.0/(pf.domain_dimensions * 2)
        xf, yf, zf = pf.domain_right_edge - 1.0/(pf.domain_dimensions * 2)
        dd = pf.h.all_data()
        rho_tot = dd.quantities["TotalQuantity"]("Density")[0]
        coords = np.mgrid[xi:xf:xn*1j, yi:yf:yn*1j, zi:zf:zn*1j]
        uc = [np.unique(c) for c in coords]
        # Some simple projection tests with single grids
        for ax, an in enumerate("xyz"):
            xax = x_dict[ax]
            yax = y_dict[ax]
            for wf in ["Density", None]:
                proj = pf.h.proj(ax, ["Ones", "Density"], weight_field = wf)
                yield assert_equal, proj["Ones"].sum(), proj["Ones"].size
                yield assert_equal, proj["Ones"].min(), 1.0
                yield assert_equal, proj["Ones"].max(), 1.0
                yield assert_equal, np.unique(proj["px"]), uc[xax]
                yield assert_equal, np.unique(proj["py"]), uc[yax]
                yield assert_equal, np.unique(proj["pdx"]), 1.0/(dims[xax]*2.0)
                yield assert_equal, np.unique(proj["pdy"]), 1.0/(dims[yax]*2.0)
            # wf == None
            yield assert_equal, wf, None
            v1 = proj["Density"].sum()
            v2 = (dd["Density"] * dd["d%s" % an]).sum()
            yield assert_rel_equal, v1, v2, 10
