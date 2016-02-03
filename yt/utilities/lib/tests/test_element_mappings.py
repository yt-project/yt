"""
This file contains tests of the intracell interpolation code contained is
yt/utilities/lib/element_mappings.pyx.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import numpy as np
from yt.testing import assert_almost_equal
from yt.utilities.lib.element_mappings import \
    test_tetra_sampler, \
    test_hex_sampler, \
    test_tri_sampler, \
    test_quad_sampler, \
    test_hex20_sampler, \
    test_wedge_sampler


def check_all_vertices(sampler, vertices, field_values):
    NV = vertices.shape[0]
    NDIM = vertices.shape[1]
    x = np.empty(NDIM)
    for i in range(NV):
        x = vertices[i]
        val = sampler(vertices, field_values, x)
        assert_almost_equal(val, field_values[i])


def test_P1Sampler2D():

    vertices = np.array([[0.1, 0.2], [0.6, 0.3], [0.2, 0.7]])
    field_values = np.array([ 1.,  2.,  3.])

    check_all_vertices(test_tri_sampler, vertices, field_values)


def test_P1Sampler3D():
    vertices = np.array([[0.1,  0.1,  0.1],
                         [0.6,  0.3,  0.2],
                         [0.2,  0.7,  0.2],
                         [0.4,  0.4,  0.7]])

    field_values = np.array([1.0, 2.0, 3.0, 4.0])

    check_all_vertices(test_tetra_sampler, vertices, field_values)


def test_Q1Sampler2D():

    vertices = np.array([[0.1, 0.2],
                         [0.6, 0.3],
                         [0.7, 0.9],
                         [0.2, 0.7]])

    field_values = np.array([1.,  2.,  3.,  4.])

    check_all_vertices(test_quad_sampler, vertices, field_values)


def test_Q1Sampler3D():
    vertices = np.array([[2.00657905, 0.6888599,  1.4375],
                         [1.8658198,  1.00973171, 1.4375],
                         [1.97881594, 1.07088163, 1.4375],
                         [2.12808879, 0.73057381, 1.4375],
                         [2.00657905, 0.6888599,  1.2   ],
                         [1.8658198,  1.00973171, 1.2   ],
                         [1.97881594, 1.07088163, 1.2   ],
                         [2.12808879, 0.73057381, 1.2   ]])

    field_values = np.array([0.4526278, 0.45262656, 0.45262657, 0.4526278,
                             0.54464296, 0.54464149, 0.5446415, 0.54464296])

    check_all_vertices(test_hex_sampler, vertices, field_values)


def test_S2Sampler3D():
    vertices = np.array([[3.00608789e-03, 4.64941000e-02, -3.95758979e-04],
                         [3.03202730e-03, 4.64941000e-02,  0.00000000e+00],
                         [3.03202730e-03, 4.70402000e-02,  3.34389809e-20],
                         [3.00608789e-03, 4.70402000e-02, -3.95758979e-04],
                         [2.45511948e-03, 4.64941000e-02, -3.23222611e-04],
                         [2.47630461e-03, 4.64941000e-02,  1.20370622e-35],
                         [2.47630461e-03, 4.70402000e-02,  3.34389809e-20],
                         [2.45511948e-03, 4.70402000e-02, -3.23222611e-04],
                         [3.01905760e-03, 4.64941000e-02, -1.97879489e-04],
                         [3.03202730e-03, 4.67671500e-02,  3.34389809e-20],
                         [3.01905760e-03, 4.70402000e-02, -1.97879489e-04],
                         [3.00608789e-03, 4.67671500e-02, -3.95758979e-04],
                         [2.73060368e-03, 4.64941000e-02, -3.59490795e-04],
                         [2.75416596e-03, 4.64941000e-02, -1.86574463e-34],
                         [2.75416596e-03, 4.70402000e-02,  6.68779617e-20],
                         [2.73060368e-03, 4.70402000e-02, -3.59490795e-04],
                         [2.47100265e-03, 4.64941000e-02, -1.61958070e-04],
                         [2.47630461e-03, 4.67671500e-02,  1.67194904e-20],
                         [2.47100265e-03, 4.70402000e-02, -1.61958070e-04],
                         [2.45511948e-03, 4.67671500e-02, -3.23222611e-04]])

    field_values = np.array([659.80151432, 650.95995348, 650.02809796,
                             658.81589888, 659.77560908, 650.93582507,
                             649.99987015, 658.78508795, 655.38073390,
                             650.49402572, 654.42199842, 659.30870660,
                             659.78856170, 650.94788928, 650.01398406,
                             658.80049342, 655.35571708, 650.46784761,
                             654.39247905, 659.28034852])

    check_all_vertices(test_hex20_sampler, vertices, field_values)


def test_W1Sampler3D():

    vertices = np.array([[-0.34641016,  0.3       ,  0.        ],
                         [-0.31754265,  0.25      ,  0.        ],
                         [-0.28867513,  0.3       ,  0.        ],
                         [-0.34641016,  0.3       ,  0.05      ],
                         [-0.31754265,  0.25      ,  0.05      ],
                         [-0.28867513,  0.3       ,  0.05      ]])

    field_values = np.array([1.,  2.,  3.,  4., 5., 6.])

    check_all_vertices(test_wedge_sampler, vertices, field_values)
