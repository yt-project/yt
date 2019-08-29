"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.utilities.answer_testing import utils


# Test data
gc = "GaussianCloud/data.0077.3d.hdf5"
tb = "TurbBoxLowRes/data.0005.3d.hdf5"
iso = "IsothermalSphere/data.0000.3d.hdf5"
zp = "ZeldovichPancake/plt32.2d.hdf5"
kho = "KelvinHelmholtz/data.0004.hdf5"


@pytest.fixture(scope='class')
def ds_gc():
    ds = utils.data_dir_load(gc)
    assert str(ds) ==  "data.0077.3d.hdf5"
    return ds

@pytest.fixture(scope='class')
def ds_tb():
    ds = utils.data_dir_load(tb)
    assert str(ds) == "data.0005.3d.hdf5"
    return ds

@pytest.fixture(scope='class')
def ds_iso():
    ds = utils.data_dir_load(iso)
    assert str(ds) == "data.0000.3d.hdf5"
    return ds

@pytest.fixture(scope='class')
def ds_zp():
    ds = utils.data_dir_load(zp)
    assert str(ds) == "plt32.2d.hdf5"
    return ds

@pytest.fixture(scope='class')
def ds_kho():
    ds = utils.data_dir_load(kho)
    assert str(ds) == "data.0004.hdf5"
    return ds
