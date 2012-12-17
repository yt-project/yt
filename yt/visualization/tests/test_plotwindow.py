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

ext_to_mime = {'.ps'  : 'application/postscript',
               '.eps' : 'application/postscript',
               '.pdf' : 'application/pdf',
               '.png' : 'image/png' }

def assert_fn(fn):
    if fn is None:
        return
    try:
        import magic
        ext = os.path.splitext(fn)[1]
        magic_text = magic.from_file(fn,mime=True)
        print fn, magic_text, ext_to_mime[ext]
        assert magic_text == ext_to_mime[ext]
    except ImportError:
        # OS X doesn't come with libmagic
        pass    

def test_plotwindow():
    pf = fake_random_pf(64)
    fns = []
    test_flnms = [None, 'test.png', 'test.eps', 
                  'test.ps', 'test.pdf']
    for fn in test_flnms:
        for dim in [0,1,2]:
            slc = SlicePlot(pf, dim, 'Density')
            fns.append(slc.save(fn)[0])
            assert_fn(fn)

            prj = ProjectionPlot(pf, dim, 'Density')
            fns.append(prj.save(fn)[0])
            assert_fn(fn)
        
        normal = [1,1,1]
        
        oaslc = OffAxisSlicePlot(pf, normal, 'Density')
        fns.append(oaslc.save(fn)[0])
        assert_fn(fn)

        oaprj = OffAxisProjectionPlot(pf, normal, 'Density')
        fns.append(oaprj.save(fn)[0])
        assert_fn(fn)
    
    teardown_func(fns)
    
