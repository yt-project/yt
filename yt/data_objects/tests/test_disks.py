from nose.tools import assert_raises
from numpy.testing import assert_equal

from yt.testing import fake_random_ds


def test_bad_disk_input():
    # Fixes 1768
    ds = fake_random_ds(16)
    with assert_raises(TypeError) as ex:
        disk = ds.disk(ds.domain_center, [0, 0, 1], ds.domain_center, (20, 'kpc'))
        disk['density']
    desired = ("Expected a numeric value (or size-1 array),"
               " received [0.5 0.5 0.5] code_length")
    assert_equal(str(ex.exception)[:50], desired[:50])

    # Valid disk
    disk = ds.disk(ds.domain_center, [0, 0, 1], (10, 'kpc'), (20, 'kpc'))
    disk['density']
