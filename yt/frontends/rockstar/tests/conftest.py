"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
r1 = "rockstar_halos/halos_0.0.bin"

_fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_mass")


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_fields_r1':
        metafunc.parametrize('f', _fields, ids=['x', 'y', 'z', 'mass'])

@pytest.fixture(scope='class')
def ds_r1():
    ds = utils.data_dir_load(r1)
    assert str(ds) == os.path.basename(r1)
    return ds
