import warnings
import numpy as np

from yt.testing import \
    assert_array_equal, \
    assert_array_almost_equal, \
    fake_random_ds
import yt.utilities.linear_interpolators as lin
from yt.utilities.lib.interpolators import \
    ghost_zone_interpolate

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
    yield assert_array_almost_equal, tfi(fv), random_data

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
    

def test_ghost_zone_extrapolation():
    ds = fake_random_ds(16)

    g = ds.index.grids[0]
    vec = g.get_vertex_centered_data(['x', 'y', 'z'], no_ghost=True)
    for i, ax in enumerate('xyz'):
        xc = g[ax]

        tf = lin.TrilinearFieldInterpolator(xc,
                (g.LeftEdge[0] + g.dds[0]/2.0,
                    g.RightEdge[0] - g.dds[0]/2.0,
                 g.LeftEdge[1] + g.dds[1]/2.0,
                    g.RightEdge[1] - g.dds[1]/2.0,
                 g.LeftEdge[2] + g.dds[2]/2.0,
                    g.RightEdge[2] - g.dds[2]/2.0),
                ["x", "y", "z"], truncate = True)

        lx,ly,lz = np.mgrid[g.LeftEdge[0]:g.RightEdge[0]:(g.ActiveDimensions[0]+1)*1j,
                            g.LeftEdge[1]:g.RightEdge[1]:(g.ActiveDimensions[1]+1)*1j,
                            g.LeftEdge[2]:g.RightEdge[2]:(g.ActiveDimensions[2]+1)*1j]
        xi  = tf({'x':lx, 'y':ly, 'z':lz})

        xz = np.zeros(g.ActiveDimensions+1)
        ghost_zone_interpolate(1, xc, np.array([0.5, 0.5, 0.5], dtype="f8"),
                                  xz, np.array([0.0, 0.0, 0.0], dtype="f8"))

        ii = (lx, ly, lz)[i]
        yield assert_array_equal, ii, vec[ax]
        yield assert_array_equal, ii, xi
        yield assert_array_equal, ii, xz


def test_get_vertex_centered_data():
    ds = fake_random_ds(16)
    g = ds.index.grids[0]

    vec_list = g.get_vertex_centered_data([('gas', 'density')], no_ghost=True)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        vec_str = g.get_vertex_centered_data('density', no_ghost=True)
        assert len(w) == 1
        assert issubclass(w[-1].category, DeprecationWarning)
        assert 'requires list of fields' in str(w[-1].message)
    vec_tuple = g.get_vertex_centered_data(('gas', 'density'), no_ghost=True) 
    assert_array_equal(vec_list[('gas', 'density')], vec_str)
    assert_array_equal(vec_list[('gas', 'density')], vec_tuple)
