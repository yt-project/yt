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
    test_tri_sampler


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
