from yt.testing import *
from yt.visualization.api import Streamlines

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

_fields = ("Density", "x-velocity", "y-velocity", "z-velocity")

def test_streamlines():
    # We decompose in different ways
    cs = np.mgrid[0.47:0.53:2j,0.47:0.53:2j,0.47:0.53:2j]
    cs = np.array([a.ravel() for a in cs]).T
    length = (1.0/128) * 16 # 16 half-widths of a cell
    for nprocs in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs = nprocs, fields = _fields)
        streams = Streamlines(pf, cs, length=length)
        streams.integrate_through_volume()
        for path in (streams.path(i) for i in range(8)):
            yield assert_rel_equal, path['dts'].sum(), 1.0, 14
            yield assert_equal, np.all(path['t'] <= (1.0 + 1e-10)), True
            path["Density"]
