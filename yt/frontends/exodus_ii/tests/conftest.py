"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    "test_displacement_fields": {
        "disp": [
            (
                {"connect2": (5.0, [0.0, 0.0, 0.0])},
                {
                    "connect1": (1.0, [1.0, 2.0, 3.0]),
                    "connect2": (0.0, [0.0, 0.0, 0.0]),
                },
            ),
            ("d1", "d2"),
        ]
    }
}


def pytest_generate_tests(metafunc):
    # Loop over each test in test_params
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            # Parametrize
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
