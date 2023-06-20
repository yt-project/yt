import numpy as np
from numpy.testing import assert_array_equal, assert_equal

import yt.utilities.lib.bitarray as ba


def test_inout_bitarray():
    # Check that we can do it for bitarrays that are funny-shaped
    for i in range(7):
        # Check we can feed in an array
        arr_in = np.random.random(32**3 + i) > 0.5
        b = ba.bitarray(arr=arr_in)
        if i > 0:
            assert_equal(b.ibuf.size, (32**3) / 8.0 + 1)
        arr_out = b.as_bool_array()
        assert_equal(arr_in, arr_out)

        # Let's check we can do it without feeding it at first
        b = ba.bitarray(size=arr_in.size)
        b.set_from_array(arr_in)
        arr_out = b.as_bool_array()
        assert_equal(arr_in, arr_out)

    # Try a big array
    arr_in = np.random.random(32**3 + i) > 0.5
    b = ba.bitarray(arr=arr_in)
    arr_out = b.as_bool_array()
    assert_equal(arr_in, arr_out)

    # Let's check we can do something interesting.
    arr_in1 = np.random.random(32**3) > 0.5
    arr_in2 = np.random.random(32**3) > 0.5
    b1 = ba.bitarray(arr=arr_in1)
    b2 = ba.bitarray(arr=arr_in2)
    b3 = ba.bitarray(arr=(arr_in1 & arr_in2))
    assert_equal((b1.ibuf & b2.ibuf), b3.ibuf)

    b = ba.bitarray(10)
    for i in range(10):
        b.set_value(i, 2)  # 2 should evaluate to True
        arr = b.as_bool_array()
        assert_equal(arr[: i + 1].all(), True)
        assert_equal(arr[i + 1 :].any(), False)
    for i in range(10):
        b.set_value(i, 0)
    arr = b.as_bool_array()
    assert_equal(arr.any(), False)
    b.set_value(7, 1)
    arr = b.as_bool_array()
    assert_array_equal(arr, [0, 0, 0, 0, 0, 0, 0, 1, 0, 0])
    b.set_value(2, 1)
    arr = b.as_bool_array()
    assert_array_equal(arr, [0, 0, 1, 0, 0, 0, 0, 1, 0, 0])


def test_set_range():
    b = ba.bitarray(127)
    # Test once where we're in the middle of start and end bits
    b.set_range(4, 65, 1)
    comparison_array = np.zeros(127, dtype="uint8")
    comparison_array[4:65] = 1
    arr = b.as_bool_array().astype("uint8")
    assert_array_equal(arr, comparison_array)

    # Test now where we're in the middle of start
    b = ba.bitarray(64)
    b.set_range(33, 36, 1)
    comparison_array = np.zeros(64, dtype="uint8")
    comparison_array[33:36] = 1
    arr = b.as_bool_array().astype("uint8")
    assert_array_equal(arr, comparison_array)

    # Now we test when we end on a byte edge, but we have 65 entries
    b = ba.bitarray(65)
    b.set_range(32, 64, 1)
    comparison_array = np.zeros(65, dtype="uint8")
    comparison_array[32:64] = 1
    arr = b.as_bool_array().astype("uint8")
    assert_array_equal(arr, comparison_array)

    # Let's do the inverse
    b = ba.bitarray(127)
    b.set_range(0, 127, 1)
    assert_equal(b.as_bool_array().all(), True)
    b.set_range(0, 127, 0)
    assert_equal(b.as_bool_array().any(), False)
    b.set_range(3, 9, 1)
    comparison_array = np.zeros(127, dtype="uint8")
    comparison_array[3:9] = 1
    arr = b.as_bool_array().astype("uint8")
    assert_array_equal(arr, comparison_array)
    b.set_range(7, 10, 0)
    comparison_array[7:10] = 0
    arr = b.as_bool_array().astype("uint8")
    assert_array_equal(arr, comparison_array)
