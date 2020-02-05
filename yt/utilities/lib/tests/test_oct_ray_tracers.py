"""Tests for octree raytracer"""

import yt

import numpy as np
from yt.utilities.lib import octree_raytracing

from itertools import product


def test_raytracing():
    ds = yt.load('ramses_deep_noamr/output_00002/info_00002.txt')

    ad = ds.all_data()
    dom = ds.index.domains[0]
    oh = dom.oct_handler

    l = 1000
    inds = octree_raytracing.domain2ind(oh, ad.selector).reshape(-1, 8)
    def get_x_octray(start, direction):
        ray = octree_raytracing.Ray(start, direction, l)

        oct_inds, cell_inds, t = octree_raytracing.ray_step(oh, ray)
        oct_inds, cell_inds, t = oct_inds[t>0], cell_inds[t>0], t[t>0]

        indexes = []
        for i in range(oct_inds.size):
            indexes.append(inds[oct_inds[i], cell_inds[i]])
        return lambda f: ad[f][indexes]

    def get_x_yt(start, direction):
        ytray = ds.ray(start, start+direction*l)

        order = ytray['t'].argsort()
        return lambda f: ytray[f][order]

    np.random.seed(16091993)
    # Test random *positive* direction.
    # NOTE: ytrays do not include the starting cell in the ray, while this implementation does

    def check(start, direction):
        oct = get_x_octray(start, direction)
        ytray = get_x_yt(start, direction)
        for field in ('x', 'y', 'z'):
            a = oct(field)
            b = ytray(field)
            np.testing.assert_allclose(a, b)
            # N = min(a.size, b.size)  # skip first few cells, see above
            # np.testing.assert_allclose(a[-N:], b[-N:])

    for i in range(100):
        start = np.random.rand(3)
        direction = np.random.rand(3)

        yield check, start, direction


def test_shallow():
    ds = yt.load('ramses_shallow/output_00002/info_00002.txt')

    ad = ds.all_data()
    dom = ds.index.domains[0]
    oh = dom.oct_handler

    l = 1000
    inds = octree_raytracing.domain2ind(oh, ad.selector).reshape(-1, 8)
    def get_x_octray(start, direction):
        ray = octree_raytracing.Ray(start, direction, l)

        oct_inds, cell_inds, t = octree_raytracing.ray_step(oh, ray)
        oct_inds, cell_inds, t = oct_inds[t>0], cell_inds[t>0], t[t>0]

        indexes = []
        for i in range(oct_inds.size):
            indexes.append(inds[oct_inds[i], cell_inds[i]])
        return lambda f: ad[f][indexes]

    def get_x_yt(start, direction):
        ytray = ds.ray(start, start+direction*l)

        order = ytray['t'].argsort()
        return lambda f: ytray[f][order]

    np.random.seed(16091993)
    # Test random *positive* direction.
    # NOTE: ytrays do not include the starting cell in the ray, while this implementation does

    def check(start, direction):
        oct = get_x_octray(start, direction)
        ytray = get_x_yt(start, direction)
        for field in ('x', 'y', 'z'):
            a = oct(field)
            b = ytray(field)
            np.testing.assert_allclose(a, b)
            # N = min(a.size, b.size)  # skip first few cells, see above
            # np.testing.assert_allclose(a[-N:], b[-N:])

    for i in range(10):
        start = np.random.rand(3)
        direction = np.random.rand(3)

        yield check, start, direction
