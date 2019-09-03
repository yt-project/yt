"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
hello_world = "hello-0210/hello-0210.block_list"
ep_cosmo = "ENZOP_DD0140/ENZOP_DD0140.block_list"


@pytest.fixture(scope='class')
def ds_hello_world():
    ds = utils.data_dir_load(hello_world)
    return ds

@pytest.fixture(scope='class')
def ds_ep_cosmo():
    ds = utils.data_dir_load(ep_cosmo)
    return ds
