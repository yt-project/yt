"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
hello_world = "hello-0210/hello-0210.block_list"
ep_cosmo = "ENZOP_DD0140/ENZOP_DD0140.block_list"


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    'test_hello_world' : {
        'f' : [('density', 'total_energy', 'velocity_x', 'velocity_y'),
                ('density', 'total_energy', 'vx', 'vy')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.25, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')]
    },
    'test_particle_fields' : {
        'f' : [('particle_position_x', 'particle_position_y', 'particle_position_z',
                'particle_velocity_x', 'particle_velocity_y', 'particle_velocity_z'),
                ('x', 'y', 'z', 'vx', 'vy', 'vz')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')]
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
def ds_hello_world():
    ds = utils.data_dir_load(hello_world)
    return ds

@pytest.fixture(scope='class')
def ds_ep_cosmo():
    ds = utils.data_dir_load(ep_cosmo)
    return ds
