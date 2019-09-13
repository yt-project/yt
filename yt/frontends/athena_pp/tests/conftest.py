"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
disk = "KeplerianDisk/disk.out1.00000.athdf"
AM06 = "AM06/AM06.out1.00400.athdf"


@pytest.fixture(scope='class')
def ds_disk():
    ds = utils.data_dir_load(disk)
    assert str(ds) == 'disk.out1.00000'
    return ds

@pytest.fixture(scope='class')
def ds_AM06():
    ds = utils.data_dir_load(AM06)
    assert str(ds) == 'AM06.out1.00400'
    return ds
