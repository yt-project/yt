"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
gc = "GaussianCloud/data.0077.3d.hdf5"
tb = "TurbBoxLowRes/data.0005.3d.hdf5"
iso = "IsothermalSphere/data.0000.3d.hdf5"
zp = "ZeldovichPancake/plt32.2d.hdf5"
kho = "KelvinHelmholtz/data.0004.hdf5"


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    'test_gc' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [('density', 'velocity_magnitude', 'magnetic_field_x'),
                ('density', 'velocity_magnitude', 'magnetic_field_x')]
    },
    'test_tb' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [('density', 'velocity_magnitude', 'magnetic_field_x'),
                ('density', 'velocity_magnitude', 'magnetic_field_x')]
    },
    'test_iso' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [('density', 'velocity_magnitude', 'magnetic_field_x'),
                ('density', 'velocity_magnitude', 'magnetic_field_x')]
    },
    'test_kho' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [('density', 'velocity_magnitude', 'magnetic_field_x'),
                ('density', 'velocity_magnitude', 'magnetic_field_x')]
    },
    'test_zp' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('c', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'rhs'), ('None', 'rhs')],
        'f' : [('rhs', 'phi'), ('rhs', 'phi')]
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
def ds_gc():
    ds = utils.data_dir_load(gc)
    assert str(ds) ==  "data.0077.3d.hdf5"
    return ds

@pytest.fixture(scope='class')
def ds_tb():
    ds = utils.data_dir_load(tb)
    assert str(ds) == "data.0005.3d.hdf5"
    return ds

@pytest.fixture(scope='class')
def ds_iso():
    ds = utils.data_dir_load(iso)
    assert str(ds) == "data.0000.3d.hdf5"
    return ds

@pytest.fixture(scope='class')
def ds_zp():
    ds = utils.data_dir_load(zp)
    assert str(ds) == "plt32.2d.hdf5"
    return ds

@pytest.fixture(scope='class')
def ds_kho():
    ds = utils.data_dir_load(kho)
    assert str(ds) == "data.0004.hdf5"
    return ds
