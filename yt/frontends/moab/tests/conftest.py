"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""


# Test data
c5 = "c5/c5.h5m"


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    "test_cantor_5": {
        "d": [
            (
                None,
                ("sphere", ("c", (0.1, "unitary"))),
                ("sphere", ("c", (0.2, "unitary"))),
            ),
            ("None", "sphere1", "sphere2"),
        ],
        "f": [(("moab", "flux"),), ("flux",)],
    }
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
