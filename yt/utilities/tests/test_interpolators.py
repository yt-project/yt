from yt.testing import *
import yt.utilities.linear_interpolators as lin

def setup():
    pass

def test_linear_interpolator_1d():
    random_data = np.random.random(64)
    fv = {'x': np.mgrid[0.0:1.0:64j]}
    # evenly spaced bins
    ufi = lin.UnilinearFieldInterpolator(random_data, (0.0, 1.0), "x", True)
    yield assert_array_equal, ufi(fv), random_data
    
    # randomly spaced bins
    size = 64
    shift = (1. / size) * np.random.random(size) - (0.5 / size)
    fv["x"] += shift
    ufi = lin.UnilinearFieldInterpolator(random_data, 
                                         np.linspace(0.0, 1.0, size) + shift, 
                                         "x", True)
    yield assert_array_almost_equal, ufi(fv), random_data, 15

def test_linear_interpolator_2d():
    random_data = np.random.random((64, 64))
    # evenly spaced bins
    fv = dict((ax, v) for ax, v in zip("xyz",
               np.mgrid[0.0:1.0:64j, 0.0:1.0:64j]))
    bfi = lin.BilinearFieldInterpolator(random_data,
            (0.0, 1.0, 0.0, 1.0), "xy", True)
    yield assert_array_equal, bfi(fv), random_data

    # randomly spaced bins
    size = 64
    bins = np.linspace(0.0, 1.0, size)
    shifts = dict((ax, (1. / size) * np.random.random(size) - (0.5 / size)) \
                  for ax in "xy")
    fv["x"] += shifts["x"][:, np.newaxis]
    fv["y"] += shifts["y"]
    bfi = lin.BilinearFieldInterpolator(random_data,
            (bins + shifts["x"], bins + shifts["y"]), "xy", True)
    yield assert_array_almost_equal, bfi(fv), random_data, 15

def test_linear_interpolator_3d():
    random_data = np.random.random((64, 64, 64))
    # evenly spaced bins
    fv = dict((ax, v) for ax, v in zip("xyz",
               np.mgrid[0.0:1.0:64j, 0.0:1.0:64j, 0.0:1.0:64j]))
    tfi = lin.TrilinearFieldInterpolator(random_data,
            (0.0, 1.0, 0.0, 1.0, 0.0, 1.0), "xyz", True)
    yield assert_array_equal, tfi(fv), random_data

    # randomly spaced bins
    size = 64
    bins = np.linspace(0.0, 1.0, size)
    shifts = dict((ax, (1. / size) * np.random.random(size) - (0.5 / size)) \
                  for ax in "xyz")
    fv["x"] += shifts["x"][:, np.newaxis, np.newaxis]
    fv["y"] += shifts["y"][:, np.newaxis]
    fv["z"] += shifts["z"]
    tfi = lin.TrilinearFieldInterpolator(random_data,
            (bins + shifts["x"], bins + shifts["y"], 
             bins + shifts["z"]), "xyz", True)
    yield assert_array_almost_equal, tfi(fv), random_data, 15
    
