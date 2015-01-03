import yt.utilities.lib.bitarray as ba
import numpy as np
from yt.testing import *

def test_inout_bitarray():
    # Check that we can do it for bitarrays that are funny-shaped
    for i in range(7):
        # Check we can feed in an array
        arr_in = (np.random.random(32**3 + i) > 0.5)
        b = ba.bitarray(arr = arr_in)
        arr_out = b.as_bool_array()
        yield assert_equal, arr_in.view("bool"), arr_out.view("bool")

        # Let's check we can do it without feeding it at first
        b = ba.bitarray(size = arr_in.size)
        b.set_from_array(arr_in)
        arr_out = b.as_bool_array()
        yield assert_equal, arr_in.view("bool"), arr_out.view("bool")

    # Let's check we can do something interesting.
    arr_in1 = (np.random.random(32**3) > 0.5)
    arr_in2 = (np.random.random(32**3) > 0.5)
    b1 = ba.bitarray(arr = arr_in1)
    b2 = ba.bitarray(arr = arr_in2)
    b3 = ba.bitarray(arr = (arr_in1 & arr_in2))
    yield assert_equal, (b1.ibuf & b2.ibuf), b3.ibuf
