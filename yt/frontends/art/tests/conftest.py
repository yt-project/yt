"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"


@pytest.fixture(scope='class')
def ds_d9p():
    ds = utils.data_dir_load(d9p)
    assert str(ds) == "10MpcBox_HartGal_csf_a0.500.d"
    return ds
