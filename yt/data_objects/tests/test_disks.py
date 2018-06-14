from nose.tools import assert_raises
from numpy.testing import assert_equal

from yt import YTQuantity
from yt.testing import fake_random_ds


def test_bad_disk_input():
    # Fixes 1768
    ds = fake_random_ds(16)

    # Test invalid 3d array
    with assert_raises(TypeError) as ex:
        disk = ds.disk(ds.domain_center, [0, 0, 1, 1], (10, 'kpc'), (20, 'kpc'))
        disk['density']
    desired = ("Expected an array of size (1,3),"
               " received <class 'list'> of length 4")
    assert_equal(str(ex.exception), desired)

    # Test invalid float
    with assert_raises(TypeError) as ex:
        disk = ds.disk(ds.domain_center, [0, 0, 1],
                       ds.domain_center, (20, 'kpc'))
        disk['density']
    desired = ("Expected a numeric value (or size-1 array),"
               " received <class 'yt.units.yt_array.YTArray'> of length 3")
    assert_equal(str(ex.exception), desired)

    # Test invalid iterable
    with assert_raises(TypeError) as ex:
        disk = ds.disk(ds.domain_center, [0, 0, 1], (10, 'kpc'),
                       (20, 'kpc'), fields=YTQuantity(1, 'kpc'))
        disk['density']
    desired = ("Expected an iterable object, received"
               " <class 'yt.units.yt_array.YTQuantity'>")
    assert_equal(str(ex.exception), desired)

    # Test invalid object
    with assert_raises(TypeError) as ex:
        disk = ds.disk(ds.domain_center, [0, 0, 1], (10, 'kpc'),
                       (20, 'kpc'), ds=ds.all_data())
        disk['density']
    desired = ("Expected an object of <class 'yt.data_objects.static_output."
               "Dataset'> type, received <class 'yt.data_objects."
               "selection_data_containers.YTRegion'>")
    assert_equal(str(ex.exception), desired)

    # Test valid disk
    disk = ds.disk(ds.domain_center, [0, 0, 1], (10, 'kpc'), (20, 'kpc'))
    disk['density']
