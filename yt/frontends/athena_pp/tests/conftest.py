"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    'test_disk' : {
        'field' : [('density', 'velocity_r'), ('density', 'velocity_r')]
    },
    'test_AM06' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [('temperature', 'density', 'velocity_magnitude', 'magnetic_field_x'),
                ('temperature', 'density', 'velocity_magnitude', 'magnetic_field_x')]
    }
}


def pytest_generate_tests(metafunc):
    # Loop over each test in test_params
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            # Parametrize
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
