from yt.testing import assert_array_equal, assert_equal
from yt.utilities.lib.allocation_container import BitmaskPool


def test_bitmask_pool():
    bmp = BitmaskPool()
    assert_equal(len(bmp), 0)
    bmp.append(100)
    assert_equal(len(bmp), 1)
    assert_equal(bmp[0].size, 100)
    bmp.append(200)
    assert_equal(len(bmp), 2)
    assert_equal(bmp[0].size, 100)
    assert_equal(bmp[1].size, 200)
    assert_equal(sum(_.size for _ in bmp.to_arrays()), 300)
    arrs = bmp.to_arrays()
    assert_equal(arrs[0].size, 100)
    assert_equal(arrs[1].size, 200)
    arrs[0][:] = 1
    arrs = bmp.to_arrays()
    assert_array_equal(arrs[0], 1)
