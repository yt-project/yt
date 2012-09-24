import numpy as np
from numpy.testing import assert_array_equal
import yt.utilities.linear_interpolators as lin

def setup():
    pass


def test_linear_interpolator():
    random_data = np.random.random(128)
    x = {"Random":np.mgrid[0.0:1.0:128j]}
    ufi = lin.UnilinearFieldInterpolator(random_data, (0.0, 1.0), "Random", True)
    assert_array_equal(ufi(x), random_data)
