import nose

from yt.config import ytcfg
from yt.testing import fake_particle_ds, requires_module

def setup_func():
    ytcfg["yt", "__withintesting"] = "True"

@requires_module('firefly')
@nose.with_setup(setup_func)
def test_firefly_JSON_object():
    ds = fake_particle_ds()
    ad = ds.all_data()
    reader = ad.create_firefly_object(
        path_to_firefly = 
        )
    reader.dumpToJSON()

