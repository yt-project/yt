"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
g1 = "owls_fof_halos/groups_001/group_001.0.hdf5"
g8 = "owls_fof_halos/groups_008/group_008.0.hdf5"


fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_mass")


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ in ['test_fields_g8', 'test_fields_g1']:
        metafunc.parametrize('field', fields, ids=['x', 'y', 'z', 'mass'])


@pytest.fixture(scope='class')
def ds_g8():
    ds = utils.data_dir_load(g8)
    assert str(ds) == os.path.basename(g8)
    return ds

@pytest.fixture(scope='class')
def ds_g1():
    ds = utils.data_dir_load(g1)
    assert str(ds) == os.path.basename(g1)
    return ds
