"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    "test_fields_g5": {
        "field": [
            (
                "particle_position_x",
                "particle_position_y",
                "particle_position_z",
                "particle_velocity_x",
                "particle_velocity_y",
                "particle_velocity_z",
                "particle_mass",
                "particle_identifier",
            ),
            ("x", "y", "z", "vx", "vy", "vz", "mass", "identifier"),
        ]
    },
    "test_fields_g42": {
        "field": [
            (
                "particle_position_x",
                "particle_position_y",
                "particle_position_z",
                "particle_velocity_x",
                "particle_velocity_y",
                "particle_velocity_z",
                "particle_mass",
                "particle_identifier",
            ),
            ("x", "y", "z", "vx", "vy", "vz", "mass", "identifier"),
        ]
    },
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
