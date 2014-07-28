import numpy as np
from yt.testing import \
    fake_random_ds, assert_equal, assert_array_less
from yt.utilities.math_utils import periodic_dist


def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

def test_point_selector():
    # generate fake amr data
    ds = fake_random_ds(64, nprocs=51)
    assert(all(ds.periodicity))

    dd = ds.all_data()
    positions = np.array([dd[ax] for ax in 'xyz'])
    delta = 0.5*np.array([dd['d'+ax] for ax in 'xyz'])
    # ensure cell centers and corners always return one and
    # only one point object
    for p in positions:
        data = ds.point(p)
        assert_equal(data["ones"].shape[0], 1)
    for p in positions - delta:
        data = ds.point(p)
        assert_equal(data["ones"].shape[0], 1)
    for p in positions + delta:
        data = ds.point(p)
        assert_equal(data["ones"].shape[0], 1)
 
def test_sphere_selector():
    # generate fake data with a number of non-cubical grids
    ds = fake_random_ds(64, nprocs=51)
    assert(all(ds.periodicity))

    # aligned tests
    spheres = [ [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.5],
                [1.0, 1.0, 1.0],
                [0.25, 0.75, 0.25] ]

    for center in spheres:
        data = ds.sphere(center, 0.25)
        # WARNING: this value has not be externally verified
        dd = ds.all_data()
        dd.set_field_parameter("center", ds.arr(center, 'code_length'))
        n_outside = (dd["radius"] >= 0.25).sum()
        assert_equal(data["radius"].size + n_outside, dd["radius"].size)

        positions = np.array([data[ax] for ax in 'xyz'])
        centers = np.tile(data.center, data['x'].shape[0]).reshape(
                          data['x'].shape[0], 3).transpose()
        dist = periodic_dist(positions, centers,
                             ds.domain_right_edge-ds.domain_left_edge,
                             ds.periodicity)
        # WARNING: this value has not been externally verified
        yield assert_array_less, dist, 0.25

def test_ellipsoid_selector():
    # generate fake data with a number of non-cubical grids
    ds = fake_random_ds(64, nprocs=51)
    assert(all(ds.periodicity))

    ellipsoids = [ [0.0, 0.0, 0.0],
                   [0.5, 0.5, 0.5],
                   [1.0, 1.0, 1.0],
                   [0.25, 0.75, 0.25] ]

    # spherical ellipsoid tests
    ratios = 3*[0.25]
    for center in ellipsoids:
        data = ds.ellipsoid(center, ratios[0], ratios[1], ratios[2], 
                              np.array([1., 0., 0.]), 0.)
        data.get_data()

        dd = ds.all_data()
        dd.set_field_parameter("center", ds.arr(center, "code_length"))
        n_outside = (dd["radius"] >= ratios[0]).sum()
        assert_equal(data["radius"].size + n_outside, dd["radius"].size)

        positions = np.array([data[ax] for ax in 'xyz'])
        centers = np.tile(data.center, 
                          data.shape[0]).reshape(data.shape[0], 3).transpose()
        dist = periodic_dist(positions, centers,
                             ds.domain_right_edge-ds.domain_left_edge,
                             ds.periodicity)
        # WARNING: this value has not been externally verified
        yield assert_array_less, dist, ratios[0]

    # aligned ellipsoid tests
    ratios = [0.25, 0.1, 0.1]
    for center in ellipsoids: 
        data = ds.ellipsoid(center, ratios[0], ratios[1], ratios[2], 
                              np.array([1., 0., 0.]), 0.)
        
        # hack to compute elliptic distance
        dist2 = np.zeros(data["ones"].shape[0])
        for i,ax in enumerate('xyz'):
            positions = np.zeros((3,data["ones"].shape[0]))
            positions[i,:] = data[ax]
            centers = np.zeros((3,data["ones"].shape[0]))
            centers[i,:] = center[i]
            dist2 += (periodic_dist(positions, centers,
                                   ds.domain_right_edge-ds.domain_left_edge,
                                   ds.periodicity)/ratios[i])**2
        # WARNING: this value has not been externally verified
        yield assert_array_less, dist2, 1.0

def test_slice_selector():
    # generate fake data with a number of non-cubical grids
    ds = fake_random_ds(64, nprocs=51)
    assert(all(ds.periodicity))

    for i,d in enumerate('xyz'):
        for coord in np.arange(0.0,1.0,0.1):
            data = ds.slice(i, coord)
            data.get_data()
            v = data[d].to_ndarray()
            yield assert_equal, data.shape[0], 64**2
            yield assert_equal, data["ones"].shape[0], 64**2
            yield assert_array_less, np.abs(v - coord), 1./128.+1e-6

def test_cutting_plane_selector():
    # generate fake data with a number of non-cubical grids
    ds = fake_random_ds(64, nprocs=51)
    assert(all(ds.periodicity))

    # test cutting plane against orthogonal plane
    for i,d in enumerate('xyz'):
        norm = np.zeros(3)
        norm[i] = 1.0

        for coord in np.arange(0, 1.0, 0.1):
            center = np.zeros(3)
            center[i] = coord

            data = ds.slice(i, coord)
            data.get_data()
            data2 = ds.cutting(norm, center)
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
