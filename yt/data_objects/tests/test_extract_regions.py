import numpy as np

from yt.testing import \
    fake_random_ds, \
    assert_equal, \
    assert_almost_equal


def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_cut_region():
    # We decompose in different ways
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(64, nprocs = nprocs,
            fields = ("density", "temperature", "velocity_x"))
        # We'll test two objects
        dd = ds.all_data()
        r = dd.cut_region( [ "obj['temperature'] > 0.5",
                             "obj['density'] < 0.75",
                             "obj['velocity_x'] > 0.25" ])
        t = ( (dd["temperature"] > 0.5 ) 
            & (dd["density"] < 0.75 )
            & (dd["velocity_x"] > 0.25 ) )
        yield assert_equal, np.all(r["temperature"] > 0.5), True
        yield assert_equal, np.all(r["density"] < 0.75), True
        yield assert_equal, np.all(r["velocity_x"] > 0.25), True
        yield assert_equal, np.sort(dd["density"][t]), np.sort(r["density"])
        yield assert_equal, np.sort(dd["x"][t]), np.sort(r["x"])
        # We are disabling these, as cutting cut regions does not presently
        # work
        #r2 = r.cut_region( [ "obj['temperature'] < 0.75" ] )
        #t2 = (r["temperature"] < 0.75)
        #yield assert_equal, np.sort(r2["temperature"]), np.sort(r["temperature"][t2])
        #yield assert_equal, np.all(r2["temperature"] < 0.75), True

        # Now we can test some projections
        dd = ds.all_data()
        cr = dd.cut_region(["obj['ones'] > 0"])
        for weight in [None, "density"]:
            p1 = ds.proj("density", 0, data_source=dd, weight_field=weight)
            p2 = ds.proj("density", 0, data_source=cr, weight_field=weight)
            for f in p1.field_data:
                yield assert_almost_equal, p1[f], p2[f]
        cr = dd.cut_region(["obj['density'] > 0.25"])
        p2 = ds.proj("density", 2, data_source=cr)
        yield assert_equal, p2["density"].max() > 0.25, True
        p2 = ds.proj("density", 2, data_source=cr, weight_field = "density")
        yield assert_equal, p2["density"].max() > 0.25, True
