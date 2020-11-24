"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    "test_datacontainer_data": {
        "field": [(("grid", "density"), ("all", "particle_mass")), ("density", "mass")]
    },
    "test_covering_grid_datacontainer_data": {
        "field": [(("grid", "density"), ("all", "particle_mass")), ("density", "mass")]
    },
    "test_arbitrary_grid_datacontainer_data": {
        "field": [(("grid", "density"), ("all", "particle_mass")), ("density", "mass")]
    },
    "test_frb_datacontainer_data": {"field": [("density",), ("density",)]},
    "test_spatial_data": {"field": [(("grid", "density"),), ("density",)]},
    "test_profile_data1": {
        "field": [("temperature", "x", "density"), ("temperature", "x", "density")]
    },
    "test_profile_data2": {
        "field": [
            ("density", "x", "temperature", "y", "cell_mass"),
            ("density", "x", "temperature", "y", "cell_mass"),
        ]
    },
    "test_nonspatial_data1": {
        "field": [
            ("region_density", "sphere_density"),
            ("region_density", "sphere_density"),
        ]
    },
    "test_nonspatial_data2": {"field": [("density",), ("density",)]},
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
