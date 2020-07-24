"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""


test_params = {
    "test_gizmo_64": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("c", (0.1, "unitary")))), ("None", "sphere")],
        "f, w": [
            (
                (("gas", "density"), None),
                (("gas", "temperature"), ("gas", "density")),
                (("gas", "metallicity"), ("gas", "denisty")),
                (("gas", "O_metallicity"), ("gas", "density")),
                (("gas", "velocity_magnitude"), None),
                (("deposit", "all_count"), None),
                (("deposit", "all_cic"), None),
                (("deposit", "PartType0_density"), None),
            ),
            (
                "density-None",
                "temperature-density",
                "metallicity-density",
                "O_metallicity-density",
                "velocity_magnitude-None",
                "count-None",
                "cic-None",
                "PartType0_density-None",
            ),
        ],
    }
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
