"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
hello_world = "hello-0210/hello-0210.block_list"
ep_cosmo = "ENZOP_DD0140/ENZOP_DD0140.block_list"


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_hello_world':
        axes = [0, 1, 2]
        dso = [ None, ("sphere", ("max", (0.25, 'unitary')))]
        weights = [None, 'density']
        _fields = ("density", "total_energy",
                   "velocity_x", "velocity_y")
        metafunc.parametrize('f', _fields, ids=['density', 'total_energy',
            'vx', 'vy'])
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', dso, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'density'])
    if metafunc.function.__name__ == 'test_particle_fields':
        dso = [ None, ("sphere", ("max", (0.1, 'unitary')))]
        _pfields = ("particle_position_x", "particle_position_y",
                    "particle_position_z", "particle_velocity_x",
                    "particle_velocity_y", "particle_velocity_z")
        metafunc.parametrize('f', _pfields, ids=['x', 'y', 'z', 'vx', 'vy', 'vz'])
        metafunc.parametrize('d', dso, ids=['None', 'sphere'])


@pytest.fixture(scope='class')
def ds_hello_world():
    ds = utils.data_dir_load(hello_world)
    return ds

@pytest.fixture(scope='class')
def ds_ep_cosmo():
    ds = utils.data_dir_load(ep_cosmo)
    return ds
