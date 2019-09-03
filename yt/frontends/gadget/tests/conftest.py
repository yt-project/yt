"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
isothermal_h5 = "IsothermalCollapse/snap_505.hdf5"

iso_kwargs = dict(bounding_box=[[-3, 3], [-3, 3], [-3, 3]])


@pytest.fixture(scope='class')
def ds_isothermal_h5():
    ds = utils.data_dir_load(isothermal_h5, kwargs=iso_kwargs)
    return ds
