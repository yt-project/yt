from yt.testing import *
from yt.mods import SlicePlot, ProjectionPlot, \
    OffAxisSlicePlot, OffAxisProjectionPlot
import os

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def teardown_func(fns):
    for fn in np.unique(fns):
        os.remove(fn)

def assert_fn(fn, fns):
    if fn is not None:
        print fn, fns[-1]
        assert fn == fns[-1]

def test_plotwindow():
    pf = fake_random_pf(64)
    fns = []
    test_flnms = [None, 'test.png', 'test.eps', 
                  'test.ps', 'test.pdf']
    for fn in test_flnms:
        for dim in [0,1,2]:
            slc = SlicePlot(pf, dim, 'Density')
            fns.append(slc.save(fn)[0])
            assert_fn(fn, fns)

            prj = ProjectionPlot(pf, dim, 'Density')
            fns.append(prj.save(fn)[0])
            assert_fn(fn, fns)
        normal = [1,1,1]
        oaslc = OffAxisSlicePlot(pf, normal, 'Density')
        fns.append(oaslc.save(fn)[0])
        assert_fn(fn,fns)

        oaprj = OffAxisProjectionPlot(pf, normal, 'Density')
        fns.append(oaprj.save(fn)[0])
        assert_fn(fn, fns)
    
    teardown_func(fns)
    
