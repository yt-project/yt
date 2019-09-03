"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
os33 = "snapshot_033/snap_033.0.hdf5"


@pytest.fixture(scope='class')
def ds_os33():
    ds = utils.data_dir_load(os33)
    return ds
