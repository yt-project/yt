from yt.testing import *
import numpy as np

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_domain_point():
    pf = fake_random_pf(16, fields = ("density"))
    p = pf.point(pf.domain_center)
