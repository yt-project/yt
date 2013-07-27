import numpy as np
from yt.testing import fake_random_pf, assert_equal, assert_array_less
from yt.utilities.math_utils import periodic_dist


def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

def test_sphere_selector():
    # generate fake data with a number of non-cubical grids
    pf = fake_random_pf(64, nprocs=51)
    assert(all(pf.periodicity))

    # aligned tests
    spheres = [ [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.5],
                [1.0, 1.0, 1.0],
                [0.25, 0.75, 0.25] ]

    for center in spheres:
        data = pf.h.sphere(center, 0.25)
        data.get_data()
        # WARNING: this value has not be externally verified
        dd = pf.h.all_data()
        dd.set_field_parameter("center", center)
        n_outside = (dd["RadiusCode"] >= 0.25).sum()
        assert_equal(data.size + n_outside, dd.size)

        positions = np.array([data[ax] for ax in 'xyz'])
        centers = np.tile(data.center, 
                          data.shape[0]).reshape(data.shape[0], 3).transpose()
        dist = periodic_dist(positions, centers,
                             pf.domain_right_edge-pf.domain_left_edge,
                             pf.periodicity)
        # WARNING: this value has not been externally verified
        yield assert_array_less, dist, 0.25

def test_ellipsoid_selector():
    # generate fake data with a number of non-cubical grids
    pf = fake_random_pf(64, nprocs=51)
    assert(all(pf.periodicity))

    ellipsoids = [ [0.0, 0.0, 0.0],
                   [0.5, 0.5, 0.5],
                   [1.0, 1.0, 1.0],
                   [0.25, 0.75, 0.25] ]

    # spherical ellipsoid tests
    ratios = 3*[0.25]
    for center in ellipsoids:
        data = pf.h.ellipsoid(center, ratios[0], ratios[1], ratios[2], 
                              np.array([1., 0., 0.]), 0.)
        data.get_data()

        dd = pf.h.all_data()
        dd.set_field_parameter("center", center)
        n_outside = (dd["RadiusCode"] >= ratios[0]).sum()
        assert_equal(data.size + n_outside, dd.size)

        positions = np.array([data[ax] for ax in 'xyz'])
        centers = np.tile(data.center, 
                          data.shape[0]).reshape(data.shape[0], 3).transpose()
        dist = periodic_dist(positions, centers,
                             pf.domain_right_edge-pf.domain_left_edge,
                             pf.periodicity)
        # WARNING: this value has not been externally verified
        yield assert_array_less, dist, ratios[0]

    # aligned ellipsoid tests
    ratios = [0.25, 0.1, 0.1]
    for center in ellipsoids: 
        data = pf.h.ellipsoid(center, ratios[0], ratios[1], ratios[2], 
                              np.array([1., 0., 0.]), 0.)
        data.get_data()
        
        # hack to compute elliptic distance
        dist2 = np.zeros(data.shape[0])
        for i,ax in enumerate('xyz'):
            positions = np.zeros((3,data.shape[0]))
            positions[i,:] = data[ax]
            centers = np.zeros((3,data.shape[0]))
            centers[i,:] = center[i]
            dist2 += (periodic_dist(positions, centers,
                                   pf.domain_right_edge-pf.domain_left_edge,
                                   pf.periodicity)/ratios[i])**2
        # WARNING: this value has not been externally verified
        yield assert_array_less, dist2, 1.0

def test_slice_selector():
    # generate fake data with a number of non-cubical grids
    pf = fake_random_pf(64, nprocs=51)
    assert(all(pf.periodicity))

    for i,d in enumerate('xyz'):
        for coord in np.arange(0,1.0,0.1):
            data = pf.h.slice(i, coord)
            data.get_data()
            assert(data.shape[0] == 64**2)
            yield assert_array_less, np.abs(data[d] - coord), 1./128.+1e-6

def test_cutting_plane_selector():
    # generate fake data with a number of non-cubical grids
    pf = fake_random_pf(64, nprocs=51)
    assert(all(pf.periodicity))

    # test cutting plane against orthogonal plane
    for i,d in enumerate('xyz'):
        norm = np.zeros(3)
        norm[i] = 1.0

        for coord in np.arange(0, 1.0, 0.1):
            center = np.zeros(3)
            center[i] = coord

            data = pf.h.slice(i, coord)
            data.get_data()
            data2 = pf.h.cutting(norm, center)
            data2.get_data()

            assert(data.shape[0] == data2.shape[0])

            cells1 = np.lexsort((data['x'],data['y'],data['z']))
            cells2 = np.lexsort((data2['x'],data2['y'],data2['z']))
            for d2 in 'xyz':
                yield assert_equal, data[d2][cells1], data2[d2][cells2]

#def test_region_selector():
#
#def test_disk_selector():
#
#def test_orthoray_selector():
#
#def test_ray_selector():
