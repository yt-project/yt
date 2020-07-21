"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    'test_snapshot_033' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('c', (0.1, 'unitary')))), ('None', 'sphere')],
        'f, w' : [((("gas", "density"), None),
                    (("gas", "temperature"), None),
                    (("gas", "temperature"), ("gas", "density")),
                    (("gas", "velocity_magnitude"), None)),
                    ("density-None", "temperature-None", "temperature-density",
                    "velocity-None")
                ]
    }
}
        

def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
