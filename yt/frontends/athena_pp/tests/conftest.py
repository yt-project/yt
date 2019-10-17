"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
disk = "KeplerianDisk/disk.out1.00000.athdf"
AM06 = "AM06/AM06.out1.00400.athdf"


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_disk':
        fields = ('density', 'velocity_r')
        metafunc.parametrize('field', fields, ids=['density', 'velocity_r'])
    if metafunc.function.__name__ == 'test_AM06':
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("temperature",
            "density",
            "velocity_magnitude",
            "magnetic_field_x"
        )
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'density'])
        metafunc.parametrize('f', fields, ids=['temperature', 'density',
            'velocity_magnitude', 'magnetic_field_x'])


@pytest.fixture(scope='class')
def ds_disk():
    ds = utils.data_dir_load(disk)
    assert str(ds) == 'disk.out1.00000'
    return ds

@pytest.fixture(scope='class')
def ds_AM06():
    ds = utils.data_dir_load(AM06)
    assert str(ds) == 'AM06.out1.00400'
    return ds
