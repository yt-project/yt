"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Data file
sizmbhloz = "sizmbhloz-clref04SNth-rs9_a0.9011/"
sizmbhloz += "sizmbhloz-clref04SNth-rs9_a0.9011.art"


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_sizmbhloz':
        dso = [None, ("sphere", ("max", (0.1, 'unitary')))]
        axes = [0, 1, 2]
        weight_fields = [None, "density"]
        fields = ("temperature", "density", "velocity_magnitude",
                   ("deposit", "all_density"), ("deposit", "all_count"))
        metafunc.parametrize('d', dso, ids=['None', 'sphere'])
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('w', weight_fields, ids=['None', 'density'])
        metafunc.parametrize('f', fields, ids=['temperature', 'density',
            'velocity_magnitude', 'dep_all_density', 'dep_all_count'])


@pytest.fixture(scope='class')
def ds_sizmbhloz():
    ds = utils.data_dir_load(sizmbhloz)
    assert str(ds) == os.path.basename(sizmbhloz)
    return ds
