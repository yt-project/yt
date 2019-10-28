"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_d9p':
        fields = (
            ("gas", "density"),
            ("gas", "temperature"),
            ("all", "particle_mass"),
            ("all", "particle_position_x")
        )
        dso = [None, ("sphere", ("max", (0.1, 'unitary')))]
        metafunc.parametrize('f', fields, ids=['density', 'temperature',
            'particle_mass', 'x'])
        metafunc.parametrize('d', dso, ids=['None', 'sphere'])
        metafunc.parametrize('a', [0, 1, 2], ids=['0', '1', '2'])
        metafunc.parametrize('w', [None, 'density'], ids=['None', 'density'])


@pytest.fixture(scope='class')
def ds_d9p():
    ds = utils.data_dir_load(d9p)
    assert str(ds) == "10MpcBox_HartGal_csf_a0.500.d"
    return ds
