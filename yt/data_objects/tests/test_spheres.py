from yt.testing import *
import numpy as np

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_domain_sphere():
    pf = fake_random_pf(16, fields = ("density"))
    sp = pf.sphere(pf.domain_center, pf.domain_width[0])
