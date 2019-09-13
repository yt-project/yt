"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
wt = "WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030"
fid_1to3_b1 = "fiducial_1to3_b1/fiducial_1to3_b1_hdf5_part_0080"
dens_turb_mag = 'DensTurbMag/DensTurbMag_hdf5_plt_cnt_0015'


@pytest.fixture(scope='class')
def ds_sloshing():
    ds = utils.data_dir_load(sloshing)
    assert str(ds) == "sloshing_low_res_hdf5_plt_cnt_0300"
    return ds

@pytest.fixture(scope='class')
def ds_wt():
    ds = utils.data_dir_load(wt)
    assert str(ds) == "windtunnel_4lev_hdf5_plt_cnt_0030"
    return ds

@pytest.fixture(scope='class')
def ds_fid_1to3_b1():
    ds = utils.data_dir_load(fid_1to3_b1)
    return ds

@pytest.fixture(scope='class')
def ds_dens_turb_mag():
    ds = utils.data_dir_load(dens_turb_mag)
    return ds
