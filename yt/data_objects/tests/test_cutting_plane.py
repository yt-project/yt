import os
import tempfile

from yt.testing import assert_equal, fake_random_ds
from yt.units.unit_object import Unit


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def teardown_func(fns):
    for fn in fns:
        try:
            os.remove(fn)
        except OSError:
            pass


def test_cutting_plane():
    fns = []
    for nprocs in [8, 1]:
        # We want to test both 1 proc and 8 procs, to make sure that
        # parallelism isn't broken
        ds = fake_random_ds(64, nprocs=nprocs)
        center = [0.5, 0.5, 0.5]
        normal = [1, 1, 1]
        cut = ds.cutting(normal, center)
        assert_equal(cut[("index", "ones")].sum(), cut[("index", "ones")].size)
        assert_equal(cut[("index", "ones")].min(), 1.0)
        assert_equal(cut[("index", "ones")].max(), 1.0)
        pw = cut.to_pw(fields=("gas", "density"))
        for p in pw.plots.values():
            tmpfd, tmpname = tempfile.mkstemp(suffix=".png")
            os.close(tmpfd)
            p.save(name=tmpname)
            fns.append(tmpname)
        for width in [(1.0, "unitary"), 1.0, ds.quan(0.5, "code_length")]:
            frb = cut.to_frb(width, 64)
            for cut_field in [("index", "ones"), ("gas", "density")]:
                fi = ds._get_field_info("unknown", cut_field)
                data = frb[cut_field]
                assert_equal(data.info["data_source"], cut.__str__())
                assert_equal(data.info["axis"], 4)
                assert_equal(data.info["field"], str(cut_field))
                assert_equal(data.units, Unit(fi.units))
                assert_equal(data.info["xlim"], frb.bounds[:2])
                assert_equal(data.info["ylim"], frb.bounds[2:])
                assert_equal(data.info["length_to_cm"], ds.length_unit.in_cgs())
                assert_equal(data.info["center"], cut.center)
    teardown_func(fns)
