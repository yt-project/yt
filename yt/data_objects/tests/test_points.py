from yt.testing import *
import numpy as np

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_domain_point():
    ds = fake_random_ds(16, fields = ("density"))
    p = ds.point(ds.domain_center)
