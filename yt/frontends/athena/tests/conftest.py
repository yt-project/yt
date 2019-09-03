"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
cloud = "ShockCloud/id0/Cloud.0050.vtk"


@pytest.fixture(scope='class')
def ds_cloud():
    ds = utils.data_dir_load(cloud)
    return ds
