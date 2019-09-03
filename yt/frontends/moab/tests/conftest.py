"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
c5 = "c5/c5.h5m"


@pytest.fixture(scope='class')
def ds_c5():
    ds = utils.data_dir_load(c5)
    assert str(ds) == "c5"
    return ds
