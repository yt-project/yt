from yt.testing import *

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_cut_region():
    # We decompose in different ways
    for nprocs in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs = nprocs,
            fields = ("Density", "Temperature", "x-velocity"))
        # We'll test two objects
        dd = pf.h.all_data()
        r = dd.cut_region( [ "obj['Temperature'] > 0.5",
                             "obj['Density'] < 0.75",
                             "obj['x-velocity'] > 0.25" ])
        t = ( (dd["Temperature"] > 0.5 ) 
            & (dd["Density"] < 0.75 )
            & (dd["x-velocity"] > 0.25 ) )
        yield assert_equal, np.all(r["Temperature"] > 0.5), True
        yield assert_equal, np.all(r["Density"] < 0.75), True
        yield assert_equal, np.all(r["x-velocity"] > 0.25), True
        yield assert_equal, np.sort(dd["Density"][t]), np.sort(r["Density"])
        yield assert_equal, np.sort(dd["x"][t]), np.sort(r["x"])
        r2 = r.cut_region( [ "obj['Temperature'] < 0.75" ] )
        t2 = (r["Temperature"] < 0.75)
        yield assert_equal, np.sort(r2["Temperature"]), np.sort(r["Temperature"][t2])
        yield assert_equal, np.all(r2["Temperature"] < 0.75), True
