"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
r1 = "rockstar_halos/halos_0.0.bin"


@pytest.fixture(scope='class')
def ds_r1():
    ds = utils.data_dir_load(r1)
    assert str(ds) == os.path.basename(r1)
    return ds
