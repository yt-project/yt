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


axes = [0, 1, 2]
center = "max"
ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
weights = [None, "density"]


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_sloshing':
        fields = ("temperature", "density", "velocity_magnitude")
    if metafunc.function.__name__ == 'test_wind_tunnel':
        fields = ("temperature", "density")
    if metafunc.function.__name__ in ['test_sloshing', 'test_wind_tunnel']:
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'density'])
        metafunc.parametrize('f', fields, ids=[str(i) for i in range(len(fields))])
    if metafunc.function.__name__ == 'test_fid_1to3_b1':
        fid_1to3_b1_fields = { 
                ("deposit", "all_density") : None,
                ("deposit", "all_count") : None,
                ("deposit", "all_cic") : None,
                ("deposit", "all_cic_velocity_x") : ("deposit", "all_cic"),
                ("deposit", "all_cic_velocity_y") : ("deposit", "all_cic"),
                ("deposit", "all_cic_velocity_z") : ("deposit", "all_cic")
        }
        metafunc.parametrize('f, w', [(k, v) for k, v in fid_1to3_b1_fields.items()],
            ids=['all_dens', 'all_count', 'all_cic', 'all_cic_vx', 'all_cic_vy', 'all_cic_vz'])


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
