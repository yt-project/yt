from yt.testing import *
import yt.utilities.linear_interpolators as lin

def setup():
    pass

def test_linear_interpolator_1d():
    random_data = np.random.random(64)
    fv = {'x': np.mgrid[0.0:1.0:64j]}
    ufi = lin.UnilinearFieldInterpolator(random_data, (0.0, 1.0), "x", True)
    assert_array_equal(ufi(fv), random_data)

def test_linear_interpolator_2d():
    random_data = np.random.random((64, 64))
    fv = dict((ax, v) for ax, v in zip("xyz",
               np.mgrid[0.0:1.0:64j, 0.0:1.0:64j]))
    bfi = lin.BilinearFieldInterpolator(random_data,
            (0.0, 1.0, 0.0, 1.0), "xy", True)
    assert_array_equal(bfi(fv), random_data)

def test_linear_interpolator_3d():
    random_data = np.random.random((64, 64, 64))
    fv = dict((ax, v) for ax, v in zip("xyz",
               np.mgrid[0.0:1.0:64j, 0.0:1.0:64j, 0.0:1.0:64j]))
    tfi = lin.TrilinearFieldInterpolator(random_data,
            (0.0, 1.0, 0.0, 1.0, 0.0, 1.0), "xyz", True)
    assert_array_equal(tfi(fv), random_data)
