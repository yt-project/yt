"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.testing import assert_equal
from yt.utilities.answer_testing import utils


# Test data
stripping = "RamPressureStripping/id0/rps.0062.vtk"


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    'test_cloud' : {
        'f' : [('scalar[0]', 'density', 'total_energy'), ('scalar0', 'density', 'total_energy')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')]
    },
    'test_blast' : {
        'f' : [('temperature', 'density', 'velocity_magnitude'),
                ('temperature', 'density', 'velocity_magnitude')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')]
    },
    'test_stripping' : {
        'f' : [('temperature', 'density', 'specific_scalar[0]'),
                ('temperature', 'density', 'specific_scalar0')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')]
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
def ds_stripping():
    uo_stripping = {"time_unit":3.086e14,
                    "length_unit":8.0236e22,
                    "mass_unit":9.999e-30*8.0236e22**3}
    ds = utils.data_dir_load(stripping, kwargs={"units_override":uo_stripping})
    assert_equal(str(ds), "rps.0062")
    return ds
