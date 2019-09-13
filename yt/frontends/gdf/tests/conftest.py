"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
sedov = "sedov/sedov_tst_0004.h5"


@pytest.fixture(scope='class')
def ds_sedov():
    ds = utils.data_dir_load(sedov)
    assert str(ds) == "sedov_tst_0004"
    return ds
