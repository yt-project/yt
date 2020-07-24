"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
ahf_halos = "ahf_halos/snap_N64L16_135.parameter"


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    "test_fields_ahf_halos": {
        "field": [
            (
                ("nbody", "particle_position_x"),
                ("nbody", "particle_position_y"),
                ("nbody", "particle_position_z"),
                ("nbody", "particle_mass"),
            ),
            ("x", "y", "z", "mass"),
        ]
    }
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])


@pytest.fixture(scope="class")
def ds_ahf_halos():
    ds = utils.data_dir_load(ahf_halos, kwargs={"hubble_constant": 0.7})
    assert str(ds) == os.path.basename(ahf_halos)
    return ds
