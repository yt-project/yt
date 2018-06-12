from nose.tools import assert_raises, assert_equal

from yt.testing import fake_random_ds
from yt.utilities.exceptions import YTInvalidArgumentType


def test_bad_disk_input():
    # Fixes 1768
    ds = fake_random_ds(16)
    with assert_raises(YTInvalidArgumentType) as ex:
        disk = ds.disk(ds.domain_center, [0, 0, 1], ds.domain_center, (20, 'kpc'))
        disk['density']
    desired = "One of the arguments to disk (selector) does not match" \
              " its signature type."
    assert_equal(str(ex.exception)[:50], desired[:50])

    # Valid disk
    disk = ds.disk(ds.domain_center, [0, 0, 1], (10, 'kpc'), (20, 'kpc'))
    disk['density']
