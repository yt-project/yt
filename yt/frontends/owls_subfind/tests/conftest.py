"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest


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
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
