from yt.testing import *
import numpy as np

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_domain_sphere():
    ds = fake_random_ds(16, fields = ("density"))
    sp = ds.sphere(ds.domain_center, ds.domain_width[0])
