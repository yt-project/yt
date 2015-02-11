import os.path
import yt
from yt.testing import \
    assert_equal, requires_file

G30 = "IsolatedGalaxy/galaxy0030/galaxy0030"

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "serialize"] = "True"

@requires_file(G30)
def test_store():
    ds = yt.load(G30)
    store = ds.parameter_filename + '.yt'
    field = "density"
    if os.path.isfile(store):
        os.remove(store)

    proj1 = ds.proj(field, "z")
    sp = ds.sphere(ds.domain_center, (4, 'kpc'))
    proj2 = ds.proj(field, "z", data_source=sp)

    proj1_c = ds.proj(field, "z")
    yield assert_equal, proj1[field], proj1_c[field]

    proj2_c = ds.proj(field, "z", data_source=sp)
    yield assert_equal, proj2[field], proj2_c[field]
