"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""


# Test data
sedov = "sedov/sedov_tst_0004.h5"


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    "test_sedov_tunnel": {
        "axis": [(0, 1, 2), ("0", "1", "2")],
        "dobj": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "weight": [(None, ("gas", "density")), ("None", "density")],
        "field": [(("gas", "density"), ("gas", "velocity_x")), ("density", "vx")],
    }
}


def pytest_generate_tests(metafunc):
    # Loop over each test in test_params
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            # Parametrize
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
