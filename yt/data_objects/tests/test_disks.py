from nose.tools import assert_raises
from numpy.testing import assert_equal

from yt import YTQuantity
from yt.testing import fake_random_ds


def test_bad_disk_input():
    # Fixes 1768
    ds = fake_random_ds(16)

    # Test invalid 3d array
    with assert_raises(TypeError) as ex:
        ds.disk(ds.domain_center, [0, 0, 1, 1], (10, "kpc"), (20, "kpc"))
    desired = "Expected an array of size (3,), received 'list' of length 4"
    assert_equal(str(ex.exception), desired)

    # Test invalid float
    with assert_raises(TypeError) as ex:
        ds.disk(ds.domain_center, [0, 0, 1], ds.domain_center, (20, "kpc"))
    desired = (
        "Expected a numeric value (or size-1 array),"
        " received 'unyt.array.unyt_array' of length 3"
    )
    assert_equal(str(ex.exception), desired)

    # Test invalid float
    with assert_raises(TypeError) as ex:
        ds.disk(ds.domain_center, [0, 0, 1], (10, 10), (20, "kpc"))
    desired = (
        "Expected a numeric value (or tuple of format (float, String)),"
        " received an inconsistent tuple '(10, 10)'."
    )
    assert_equal(str(ex.exception), desired)

    # Test invalid iterable
    with assert_raises(TypeError) as ex:
        ds.disk(
            ds.domain_center,
            [0, 0, 1],
            (10, "kpc"),
            (20, "kpc"),
            fields=YTQuantity(1, "kpc"),
        )
    desired = "Expected an iterable object, received 'unyt.array.unyt_quantity'"
    assert_equal(str(ex.exception), desired)

    # Test invalid object
    with assert_raises(TypeError) as ex:
        ds.disk(ds.domain_center, [0, 0, 1], (10, "kpc"), (20, "kpc"), ds=ds.all_data())
    desired = (
        "Expected an object of 'yt.data_objects.static_output.Dataset' "
        "type, received "
        "'yt.data_objects.selection_objects.region.YTRegion'"
    )
    assert_equal(str(ex.exception), desired)

    # Test valid disk
    ds.disk(ds.domain_center, [0, 0, 1], (10, "kpc"), (20, "kpc"))
    ds.disk(ds.domain_center, [0, 0, 1], 10, (20, "kpc"))
