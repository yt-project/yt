"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
sedov = "sedov/sedov_tst_0004.h5"


axes = [0, 1, 2]
center = "max"
ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
weights = [None, "density"]
fields = ("density", "velocity_x")


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_sedov_tunnel':
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'density'])
        metafunc.parametrize('f', fields, ids=['density', 'vx'])


@pytest.fixture(scope='class')
def ds_sedov():
    ds = utils.data_dir_load(sedov)
    assert str(ds) == "sedov_tst_0004"
    return ds
