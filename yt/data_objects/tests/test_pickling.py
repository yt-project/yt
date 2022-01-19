import pickle

from yt.testing import assert_equal, requires_file
from yt.utilities.answer_testing.framework import data_dir_load

tipsy_ds = "TipsyGalaxy/galaxy.00300"
enzo_ds = "enzo_tiny_cosmology/DD0000/DD0000"


@requires_file(enzo_ds)
def test_grid_pickles():
    ds = data_dir_load(enzo_ds)
    ad = ds.all_data()
    # just test ad since there is a nested ds that will get (un)pickled
    ad_pickle = pickle.loads(pickle.dumps(ad))
    assert_equal(ad.ds.basename, ad_pickle.ds.basename)


@requires_file(tipsy_ds)
def test_particle_pickles():
    ds = data_dir_load(tipsy_ds)
    ad = ds.all_data()
    ds.index._identify_base_chunk(ad)
    ch0 = list(ds.index._chunk_io(ad, cache=False))[0]
    # just test one chunk since there is a nested ds that will get (un)pickled
    ch_pickle = pickle.loads(pickle.dumps(ch0))
    assert_equal(ch0.dobj.ds.basename, ch_pickle.dobj.ds.basename)
