import pytest

from yt.testing import _geom_transforms

geom_ids = [g.__repr__ for g in sorted(_geom_transforms)]


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_axial_pixelization':
        metafunc.parametrize('geom', sorted(_geom_transforms), ids=geom_ids)
        metafunc.parametrize('axis', ['x_axis', 'y_axis'], ids=['x', 'y'])
