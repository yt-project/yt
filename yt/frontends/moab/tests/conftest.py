"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
c5 = "c5/c5.h5m"

dso = [ None, ("sphere", ("c", (0.1, 'unitary'))),
          ("sphere", ("c", (0.2, 'unitary')))]
_fields = (("moab", "flux"),)


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_cantor_5':
        metafunc.parametrize('f', _fields, ids=['flux'])
        metafunc.parametrize('d', dso, ids=['None', 'sphere1', 'sphere2'])

@pytest.fixture(scope='class')
def ds_c5():
    ds = utils.data_dir_load(c5)
    assert str(ds) == "c5"
    return ds
