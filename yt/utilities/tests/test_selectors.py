from yt.testing import *
from yt.utilities.math_utils import periodic_dist

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_sphere_selector():
    # generate fake data with a number of non-cubical grids
    pf = fake_random_pf(64,nprocs=51)
    assert(all(pf.periodicity))

    # aligned tests
    spheres = [ [0.0,0.0,0.0],
                [0.5,0.5,0.5],
                [1.0,1.0,1.0],
                [0.25,0.75,0.25] ]

    for center in spheres :
        data = pf.h.sphere(center, 0.25)
        data.get_data()
        # WARNING: this value has not be externally verified
        yield assert_equal, data.size, 19568

        positions = np.array([data[ax] for ax in 'xyz'])
        centers = np.tile( data.center, data.shape[0] ).reshape(data.shape[0],3).transpose()
        dist = periodic_dist(positions, centers,
                         pf.domain_right_edge-pf.domain_left_edge,
                         pf.periodicity)
        # WARNING: this value has not been externally verified
        yield assert_almost_equal, dist.max(), 0.261806188752
