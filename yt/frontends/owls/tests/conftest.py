"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
os33 = "snapshot_033/snap_033.0.hdf5"

# This maps from field names to weight field names to use for projections
_fields = [
        [("gas", "density"), None],
        [("gas", "temperature"), None],
        [("gas", "temperature"), ("gas", "density")],
        [('gas', 'He_p0_number_density'), None],
        [('gas', 'velocity_magnitude'), None],
        [("deposit", "all_density"), None],
        [("deposit", "all_count"), None],
        [("deposit", "all_cic"), None],
        [("deposit", "PartType0_density"), None],
        [("deposit", "PartType4_density"), None]
    ]
dso = [None, ('sphere', ('c', (0.1, 'unitary')))]


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_snapshot_033':
        metafunc.parametrize('f,w', [(pair[0], pair[1]) for pair in _fields],
            ids=['dens', 'temp-None', 'temp-dens', 'He_p0_number_density',
            'velocity_magnitude', 'all_density', 'all_cic', 'PartType0Dens',
            'PartType4Dens'])
        metafunc.parametrize('d', dso, ids=['None', 'sphere'])
        metafunc.parametrize('a', [0, 1, 2], ids=['0', '1', '2'])


@pytest.fixture(scope='class')
def ds_os33():
    ds = utils.data_dir_load(os33)
    return ds
