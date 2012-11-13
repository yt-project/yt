from yt.testing import *

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_slice():
    for nprocs in [8, 1]:
        # We want to test both 1 proc and 8 procs, to make sure that
        # parallelism isn't broken
        pf = fake_random_pf(64, nprocs = nprocs)
        dims = pf.domain_dimensions
        xn, yn, zn = pf.domain_dimensions
        xi, yi, zi = pf.domain_left_edge + 1.0/(pf.domain_dimensions * 2)
        xf, yf, zf = pf.domain_right_edge - 1.0/(pf.domain_dimensions * 2)
        coords = np.mgrid[xi:xf:xn*1j, yi:yf:yn*1j, zi:zf:zn*1j]
        uc = [np.unique(c) for c in coords]
        slc_pos = 0.5
        # Some simple slice tests with single grids
        for ax, an in enumerate("xyz"):
            xax = x_dict[ax]
            yax = y_dict[ax]
            for wf in ["Density", None]:
                slc = pf.h.slice(ax, slc_pos, ["Ones", "Density"])
                yield assert_equal, slc["Ones"].sum(), slc["Ones"].size
                yield assert_equal, slc["Ones"].min(), 1.0
                yield assert_equal, slc["Ones"].max(), 1.0
                yield assert_equal, np.unique(slc["px"]), uc[xax]
                yield assert_equal, np.unique(slc["py"]), uc[yax]
                yield assert_equal, np.unique(slc["pdx"]), 1.0/(dims[xax]*2.0)
                yield assert_equal, np.unique(slc["pdy"]), 1.0/(dims[yax]*2.0)
                frb = slc.to_frb((1.0,'unitary'), 64)
                for slc_field in ['Ones', 'Density']:
                    yield assert_equal, frb[slc_field].info['data_source'], \
                            slc.__str__()
                    yield assert_equal, frb[slc_field].info['axis'], \
                            ax
                    yield assert_equal, frb[slc_field].info['field'], \
                            slc_field
                    yield assert_equal, frb[slc_field].info['units'], \
                            pf.field_info[slc_field].get_units()
                    yield assert_equal, frb[slc_field].info['xlim'], \
                            frb.bounds[:2]
                    yield assert_equal, frb[slc_field].info['ylim'], \
                            frb.bounds[2:]
                    yield assert_equal, frb[slc_field].info['length_to_cm'], \
                            pf['cm']
                    yield assert_equal, frb[slc_field].info['center'], \
                            slc.center
                    yield assert_equal, frb[slc_field].info['coord'], \
                            slc_pos
            # wf == None
            yield assert_equal, wf, None


