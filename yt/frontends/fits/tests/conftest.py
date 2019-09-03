"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import os

import pytest

from yt.frontends.fits.data_structures import FITSDataset, \
    SpectralCubeFITSDataset, \
    SkyDataFITSDataset, \
    EventsFITSDataset
from yt.utilities.answer_testing import utils


# Test data
grs = "radio_fits/grs-50-cube.fits"
vf = "UnigridData/velocity_field_20.fits"
acis = "xray_fits/acisf05356N003_evt2.fits"
A2052 = "xray_fits/A2052_merged_0.3-2_match-core_tmap_bgecorr.fits"


@pytest.fixture(scope='class')
def ds_grs():
    ds = utils.data_dir_load(grs, cls=SpectralCubeFITSDataset, kwargs={"nan_mask":0.0})
    assert str(ds) == "grs-50-cube.fits"
    return ds

@pytest.fixture(scope='class')
def ds_vf():
    ds = utils.data_dir_load(vf, cls=FITSDataset)
    assert str(ds) == "velocity_field_20.fits"
    return ds

@pytest.fixture(scope='class')
def ds_acis():
    ds = utils.data_dir_load(acis, cls=EventsFITSDataset)
    return ds

@pytest.fixture(scope='class')
def ds_A2052():
    ds = utils.data_dir_load(A2052, cls=SkyDataFITSDataset)
    assert str(ds) == "A2052_merged_0.3-2_match-core_tmap_bgecorr.fits"
    return ds
