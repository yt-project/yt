"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Data file
sizmbhloz = "sizmbhloz-clref04SNth-rs9_a0.9011/"
sizmbhloz += "sizmbhloz-clref04SNth-rs9_a0.9011.art"


@pytest.fixture(scope='class')
def ds_sizmbhloz():
    ds = utils.data_dir_load(sizmbhloz)
    assert str(ds) == os.path.basename(sizmbhloz)
    return ds
