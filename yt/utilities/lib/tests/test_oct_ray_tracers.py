"""Tests for octree raytracer"""

import yt

import numpy as np
from yt.utilities.lib import octree_raytracing


def test_raytracing():
    ds = yt.load('ramses_deep_noamr/output_00002/info_00002.txt')

    ad = ds.all_data()
    dom = ds.index.domains[0]
    oh = dom.oct_handler

    start = np.array([0.01, 1e-3, 1e-3], float)
    direction = np.array([1, 0, 0], float)
    l = 10
    direction /= np.linalg.norm(direction)

    l = 1000
    def get_x_octray(start, direction):
        ray = octree_raytracing.Ray(start, direction, l)

        oct_inds, cell_inds = octree_raytracing.ray_step(oh, ray, 1)
        inds = octree_raytracing.domain2ind(oh, ad.selector).reshape(-1, 8)

        indexes = []
        for i in range(oct_inds.size):
            indexes.append(inds[oct_inds[i], cell_inds[i]])
        return lambda f: ad[f][indexes]

    def get_x_yt(start, direction):
        ytray = ds.ray(start, start+direction*l)

        order = ytray['t'].argsort()
        return lambda f: ytray[f][order]

    np.random.seed(16091992)
    # Test random *positive* direction.
    # NOTE: ytrays do not include the starting cell in the ray, while this implementation does 
    for i in range(10):
        direction = np.random.rand(3)
        start = np.random.rand(3)
        #start[np.random.randint(3)] = 0
        oct = get_x_octray(start, direction)
        ytray = get_x_yt(start, direction)
        for field in ('x', 'y', 'z'):
            a = oct(field)
            b = ytray(field)
            N = min(a.size, b.size)  # skip first few cells, see above
            np.testing.assert_allclose(a[-N:], b[-N:])