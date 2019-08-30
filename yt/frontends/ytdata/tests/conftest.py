"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"


@pytest.fixture(scope='class')
def ds_enzotiny():
    ds = utils.data_dir_load(enzotiny)
    return ds
