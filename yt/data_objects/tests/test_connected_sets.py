from yt.utilities.answer_testing.level_sets_tests import \
    ExtractConnectedSetsTest
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load

g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"
@requires_ds(g30, big_data=True)
def test_connected_sets():
    ds = data_dir_load(g30)
    data_source = ds.disk([0.5, 0.5, 0.5], [0., 0., 1.],
                          (8, 'kpc'), (1, 'kpc'))
    yield ExtractConnectedSetsTest(g30, data_source, ("gas", "density"),
                                   5, 1e-24, 8e-24)
