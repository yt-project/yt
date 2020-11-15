from yt.testing import assert_equal, assert_raises, fake_particle_ds
from yt.utilities.exceptions import YTAmbiguousFieldName


def test_ambiguous_fails():
    ds = fake_particle_ds()
    with assert_raises(YTAmbiguousFieldName) as ex:
        _ = ds.r[:]["particle_position_x"]
    assert_equal(ex.exception.fname, "particle_position_x")
    assert_equal(ex.exception.possible_ftypes, set(("all", "nbody")))
