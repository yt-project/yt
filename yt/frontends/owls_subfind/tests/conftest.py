"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
g1 = "owls_fof_halos/groups_001/group_001.0.hdf5"
g8 = "owls_fof_halos/groups_008/group_008.0.hdf5"


@pytest.fixture(scope='class')
def ds_g8():
    ds = utils.data_dir_load(g8)
    assert str(ds) == os.path.basename(g8)
    return ds

@pytest.fixture(scope='class')
def ds_g1():
    ds = utils.data_dir_load(g1)
    assert str(ds) == os.path.basename(g1)
    return ds
