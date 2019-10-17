"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
out = "ExodusII/out.e"
out_s002 = "ExodusII/out.e-s002"
gold = "ExodusII/gold.e"
big_data = "MOOSE_sample_data/mps_out.e"


def pytest_generate_test(metafunc):
    if metafunc.function.__name__ == 'test_displacement_fields':
        displacement_dicts =[{'connect2': (5.0, [0.0, 0.0, 0.0])},
                             {'connect1': (1.0, [1.0, 2.0, 3.0]),
                              'connect2': (0.0, [0.0, 0.0, 0.0])}]
        metafunc.parametrize('disp', displacement_dicts, ids=['d1', 'd2'])


@pytest.fixture(scope='class')
def ds_out():
    ds = utils.data_dir_load(out)
    assert str(ds) == 'out.e'
    return ds

@pytest.fixture(scope='class')
def ds_out_s002():
    ds = utils.data_dir_load(out_s002)
    assert str(ds) ==  "out.e-s002"
    return ds

@pytest.fixture(scope='class')
def ds_gold():
    ds = utils.data_dir_load(gold)
    assert str(ds) == "gold.e"
    return ds
