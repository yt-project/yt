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


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    'test_fields_g8' : {
        'field' : [('particle_position_x', 'particle_position_y', 'particle_position_z',
                    'particle_mass'), ('x', 'y', 'z', 'mass')]
    },
    'test_fields_g1' : {
        'field' : [('particle_position_x', 'particle_position_y', 'particle_position_z',
                    'particle_mass'), ('x', 'y', 'z', 'mass')]
    }
}


def pytest_generate_tests(metafunc):
    # Loop over each test in test_params
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            # Parametrize
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])


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
