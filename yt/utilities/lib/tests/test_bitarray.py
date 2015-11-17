import numpy as np

import yt.utilities.lib.bitarray as ba
from yt.testing import assert_equal

def test_inout_bitarray():
    # Check that we can do it for bitarrays that are funny-shaped
    for i in range(7):
        # Check we can feed in an array
        arr_in = (np.random.random(32**3 + i) > 0.5)
        b = ba.bitarray(arr = arr_in)
        if i > 0:
            yield assert_equal, b.ibuf.size, (32**3)/8.0 + 1
        arr_out = b.as_bool_array()
        yield assert_equal, arr_in, arr_out

        # Let's check we can do it without feeding it at first
        b = ba.bitarray(size = arr_in.size)
        b.set_from_array(arr_in)
        arr_out = b.as_bool_array()
        yield assert_equal, arr_in, arr_out

    # Try a big array
    arr_in = (np.random.random(32**3 + i) > 0.5)
    b = ba.bitarray(arr = arr_in)
    arr_out = b.as_bool_array()
    yield assert_equal, arr_in, arr_out

    # Let's check we can do something interesting.
    arr_in1 = (np.random.random(32**3) > 0.5)
    arr_in2 = (np.random.random(32**3) > 0.5)
    b1 = ba.bitarray(arr = arr_in1)
    b2 = ba.bitarray(arr = arr_in2)
    b3 = ba.bitarray(arr = (arr_in1 & arr_in2))
    yield assert_equal, (b1.ibuf & b2.ibuf), b3.ibuf

    b = ba.bitarray(10)
    for i in range(10):
        b.set_value(i, 2) # 2 should evaluate to True
        arr = b.as_bool_array()
        yield assert_equal, arr[:i+1].all(), True
        yield assert_equal, arr[i+1:].any(), False
