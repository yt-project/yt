"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
from yt.utilities.answer_testing import utils

toro1d = utils.get_parameterization("ToroShockTube/DD0001/data0001")
kh2d = utils.get_parameterization("EnzoKelvinHelmholtz/DD0011/DD0011")


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    "test_toro1d": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": toro1d,
    },
    "test_kh2d": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": kh2d,
    },
    "test_moving7": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": [
            ("temperature", "density", "velocity_magnitude", "velocity_divergence"),
            ("temperature", "density", "velocity_magnitude", "velocity_divergence"),
        ],
    },
    "test_galaxy0030": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": [
            ("temperature", "density", "velocity_magnitude", "velocity_divergence"),
            ("temperature", "density", "velocity_magnitude", "velocity_divergence"),
        ],
    },
    "test_simulated_halo_mass_function": {"finder": [("fof", "hop"), ("fof", "hop")]},
    "test_analytic_halo_mass_function": {
        "fit": [[i for i in range(1, 6)], [str(i) for i in range(1, 6)]]
    },
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
