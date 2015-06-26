import numpy as np

from yt.testing import *
from yt.utilities.lib.element_mappings import \
    P1Sampler2D, \
    P1Sampler3D, \
    Q1Sampler2D, \
    Q1Sampler3D


def setup():
    pass


def test_P1Sampler2D():
    NV = 3
    NDIM = 2
    vertices = np.empty((NV, NDIM))

    vertices[0, 0] = 0.1
    vertices[0, 1] = 0.2

    vertices[1, 0] = 0.6
    vertices[1, 1] = 0.3

    vertices[2, 0] = 0.2
    vertices[2, 1] = 0.7

    field_values = np.empty(NV)
    field_values[0] = 1.0
    field_values[1] = 2.0
    field_values[2] = 3.0

    physical_x = np.empty(NDIM)
    for i in range(NV):
        physical_x = vertices[i]
        sampler = P1Sampler2D()
        x = sampler.map_real_to_unit(physical_x, vertices)
        val = P1Sampler2D.sample_at_unit_point(x, field_values)
        assert_almost_equal(field_values[i], val)


def test_P1Sampler3D():
    NV = 4
    NDIM = 3
    vertices = np.empty((NV, NDIM))

    vertices[0, 0] = 0.1
    vertices[0, 1] = 0.1
    vertices[0, 2] = 0.1

    vertices[1, 0] = 0.6
    vertices[1, 1] = 0.3
    vertices[1, 2] = 0.2

    vertices[2, 0] = 0.2
    vertices[2, 1] = 0.7
    vertices[2, 2] = 0.2

    vertices[3, 0] = 0.4
    vertices[3, 1] = 0.4
    vertices[3, 2] = 0.7

    field_values = np.empty(NV)
    field_values[0] = 1.0
    field_values[1] = 2.0
    field_values[2] = 3.0
    field_values[3] = 4.0

    physical_x = np.empty(NDIM)
    for i in range(NV):
        physical_x = vertices[i]
        sampler = P1Sampler3D()
        x = sampler.map_real_to_unit(physical_x, vertices)
        val = P1Sampler3D.sample_at_unit_point(x, field_values)
        assert_almost_equal(field_values[i], val)


def test_Q1Sampler2D():
    NV = 4
    NDIM = 2
    vertices = np.empty((NV, NDIM))

    vertices[0, 0] = 0.1
    vertices[0, 1] = 0.2

    vertices[1, 0] = 0.6
    vertices[1, 1] = 0.3

    vertices[2, 0] = 0.2
    vertices[2, 1] = 0.7

    vertices[3, 0] = 0.7
    vertices[3, 1] = 0.9

    field_values = np.empty(NV)
    field_values[0] = 1.0
    field_values[1] = 2.0
    field_values[2] = 3.0
    field_values[3] = 4.0

    physical_x = np.empty(NDIM)
    for i in range(NV):
        physical_x = vertices[i]
        sampler = Q1Sampler2D()
        x = sampler.map_real_to_unit(physical_x, vertices)
        val = Q1Sampler2D.sample_at_unit_point(x, field_values)
        assert_almost_equal(field_values[i], val)


def test_Q1Sampler3D():
    NV = 8
    NDIM = 3

    vertices = np.array([[2.00657905, 0.6888599,  1.4375],
                         [1.8658198,  1.00973171, 1.4375],
                         [1.97881594, 1.07088163, 1.4375],
                         [2.12808879, 0.73057381, 1.4375],
                         [2.00657905, 0.6888599,  1.2],
                         [1.8658198,  1.00973171, 1.2],
                         [1.97881594, 1.07088163, 1.2],
                         [2.12808879, 0.73057381, 1.2]])

    field_values = np.array([0.4526278, 0.45262656, 0.45262657, 0.4526278,
                             0.54464296, 0.54464149, 0.5446415, 0.54464296])

    physical_x = np.empty(NDIM)
    for i in range(NV):
        physical_x = vertices[i]
        sampler = Q1Sampler3D()
        x = sampler.map_real_to_unit(physical_x, vertices)
        val = Q1Sampler3D.sample_at_unit_point(x, field_values)
        assert_almost_equal(field_values[i], val)
