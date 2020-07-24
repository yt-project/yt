"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""


fields = {
    ("all", "particle_mass"): None,
    ("all", "particle_ones"): None,
    ("all", "particle_velocity_x"): ("all", "particle_mass"),
    ("all", "particle_velocity_y"): ("all", "particle_mass"),
    ("all", "particle_velocity_z"): ("all", "particle_mass"),
}

# Is a list instead of a dictionary because ("gas", "temperature")
# appeared as a duplicate key
tg_sph_fields = [
    [("gas", "density"), None],
    [("gas", "temperature"), None],
    [("gas", "temperature"), ("gas", "density")],
    [("gas", "velocity_magnitude"), None],
    [("gas", "Fe_fraction"), None],
]

tg_nbody_fields = {
    ("Stars", "Metals"): None,
}

# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    "test_pkdgrav": {
        "f, w": [
            [(k, v) for k, v in fields.items()],
            ("mass-None", "ones-None", "vx-mass", "vy-mass", "vz-mass"),
        ],
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("c", (0.3, "unitary")))), ("None", "sphere")],
    },
    "test_gasoline_dmonly": {
        "f, w": [
            [(k, v) for k, v in fields.items()],
            ("mass-None", "ones-None", "vx-mass", "vy-mass", "vz-mass"),
        ],
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("c", (0.3, "unitary")))), ("None", "sphere")],
    },
    "test_tipsy_galaxy_sph": {
        "f, w": [
            [(i[0], i[1]) for i in tg_sph_fields],
            ("dens-None", "temp-None", "temp-dens", "v-None", "Fe-None"),
        ],
        "d": [(None, ("sphere", ("c", (0.1, "unitary")))), ("None", "sphere")],
        "a": [(0, 1, 2), ("0", "1", "2")],
    },
    "test_tipsy_galaxy_nbody": {
        "f, w": [[(k, v) for k, v in tg_nbody_fields.items()], ("metals-None",)],
        "d": [(None, ("sphere", ("c", (0.1, "unitary")))), ("None", "sphere")],
        "a": [(0, 1, 2), ("0", "1", "2")],
    },
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
