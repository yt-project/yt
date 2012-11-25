from yt.testing import *
from yt.mods import SlicePlot, ProjectionPlot, \
    OffAxisSlicePlot, OffAxisProjectionPlot
import glob, os

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def teardown_func(fns):
    for fn in fns:
        os.remove(fn)

def test_plotwindow():
    pf = fake_random_pf(64)
    fns = []
    for dim in [0,1,2]:
        slc = SlicePlot(pf, dim, 'Density')
        fns.append(slc.save()[0])
        prj = ProjectionPlot(pf, dim, 'Density')
        fns.append(prj.save()[0])
    normal = [1,1,1]
    oaslc = OffAxisSlicePlot(pf, normal, 'Density')
    fns.append(oaslc.save()[0])
    oaprj = OffAxisProjectionPlot(pf, normal, 'Density')
    fns.append(oaprj.save()[0])
    teardown_func(fns)
    
