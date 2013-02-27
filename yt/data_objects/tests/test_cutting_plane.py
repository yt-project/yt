from yt.testing import *
import os

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def teardown_func(fns):
    for fn in fns:
        os.remove(fn)

def test_cutting_plane():
    for nprocs in [8, 1]:
        # We want to test both 1 proc and 8 procs, to make sure that
        # parallelism isn't broken
        pf = fake_random_pf(64, nprocs = nprocs)
        dims = pf.domain_dimensions
        center = [0.5,0.5,0.5]
        normal = [1,1,1]
        fns = []
        cut = pf.h.cutting(normal, center)
        yield assert_equal, cut["Ones"].sum(), cut["Ones"].size
        yield assert_equal, cut["Ones"].min(), 1.0
        yield assert_equal, cut["Ones"].max(), 1.0
        pw = cut.to_pw()
        fns += pw.save()
        frb = cut.to_frb((1.0,'unitary'), 64)
        for cut_field in ['Ones', 'Density']:
            yield assert_equal, frb[cut_field].info['data_source'], \
                cut.__str__()
            yield assert_equal, frb[cut_field].info['axis'], \
                4
            yield assert_equal, frb[cut_field].info['field'], \
                cut_field
            yield assert_equal, frb[cut_field].info['units'], \
                pf.field_info[cut_field].get_units()
            yield assert_equal, frb[cut_field].info['xlim'], \
                frb.bounds[:2]
            yield assert_equal, frb[cut_field].info['ylim'], \
                frb.bounds[2:]
            yield assert_equal, frb[cut_field].info['length_to_cm'], \
                pf['cm']
            yield assert_equal, frb[cut_field].info['center'], \
                cut.center
        teardown_func(fns)
