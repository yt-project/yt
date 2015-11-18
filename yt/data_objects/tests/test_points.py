from yt.testing import fake_random_ds

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_domain_point():
    ds = fake_random_ds(16)
    p = ds.point(ds.domain_center)
    p['density']
