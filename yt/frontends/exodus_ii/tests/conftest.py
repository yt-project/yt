"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
out = "ExodusII/out.e"
out_s002 = "ExodusII/out.e-s002"
gold = "ExodusII/gold.e"
big_data = "MOOSE_sample_data/mps_out.e"


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
