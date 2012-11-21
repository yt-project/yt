from yt.testing import *
from yt.mods import SlicePlot, ProjectionPlot, \
    OffAxisSlicePlot, OffAxisProjectionPlot
import glob, os

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def teardown_func():
    fns = glob.glob("*.png")
    for fn in fns:
        os.remove(fn)

def test_plotwindow():
    pf = fake_random_pf(64)
    for dim in [0,1,2]:
        slc = SlicePlot(pf, dim, 'Density')
        prj = ProjectionPlot(pf, dim, 'Density')
        slc.save()
        prj.save()
    normal = [1,1,1]
    oaslc = OffAxisSlicePlot(pf, normal, 'Density')
    oaprj = OffAxisProjectionPlot(pf, normal, 'Density')
    teardown_func()
    
