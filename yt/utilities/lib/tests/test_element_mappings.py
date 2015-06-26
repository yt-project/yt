import numpy as np

from yt.testing import *
from yt.utilities.lib.element_mappings import \
    P1Mapping2D, \
    P1Mapping3D, \
    Q1Mapping2D, \
    Q1Mapping3D


def setup():
    pass


def test_P1Mapping2D():
    NV = 3
    NDIM = 2
    vertices = np.empty((NV, NDIM))

    vertices[0, 0] = 0.1
    vertices[0, 1] = 0.2

    vertices[1, 0] = 0.6
    vertices[1, 1] = 0.3

    vertices[2, 0] = 0.2
    vertices[2, 1] = 0.7

    physical_x[0] = 0.4
    physical_x[1] = 0.4

    field_values = np.empty(NV)
    field_values[0] = 1.0
    field_values[1] = 2.0
    field_values[2] = 3.0

    physical_x = np.empty(NDIM)
    for i in range(NV):
        physical_x = vertices[i]


def test_P1Mapping3D():
    pass


def test_Q1Mapping2D():
    pass


def test_Q1Mapping3D():
    pass
