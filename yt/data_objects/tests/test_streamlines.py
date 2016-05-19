import numpy as np

from yt.testing import \
    fake_random_ds, \
    assert_rel_equal, \
    assert_equal
from yt.visualization.streamlines import Streamlines

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

_fields = ("density", "velocity_x", "velocity_y", "velocity_z")

def test_covering_grid():
    return
    # We decompose in different ways
    cs = np.mgrid[0.47:0.53:2j,0.47:0.53:2j,0.47:0.53:2j]
    cs = np.array([a.ravel() for a in cs]).T
    length = (1.0/128) * 16 # 16 half-widths of a cell
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(64, nprocs = nprocs, fields = _fields)
        streams = Streamlines(ds, cs, length=length)
        streams.integrate_through_volume()
        for path in (streams.path(i) for i in range(8)):
            yield assert_rel_equal, path['dts'].sum(), 1.0, 14
            yield assert_equal, np.all(path['t'] <= (1.0 + 1e-10)), True
            path["density"]
