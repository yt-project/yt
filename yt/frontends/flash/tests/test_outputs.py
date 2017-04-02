"""
FLASH frontend tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.testing import \
    assert_equal, \
    requires_file, \
    units_override_check
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    small_patch_amr, \
    data_dir_load, \
    sph_answer
from yt.frontends.flash.api import FLASHDataset, \
    FLASHParticleDataset
from collections import OrderedDict

_fields = ("temperature", "density", "velocity_magnitude")

sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
@requires_ds(sloshing, big_data=True)
def test_sloshing():
    ds = data_dir_load(sloshing)
    yield assert_equal, str(ds), "sloshing_low_res_hdf5_plt_cnt_0300"
    for test in small_patch_amr(ds, _fields):
        test_sloshing.__name__ = test.description
        yield test

_fields_2d = ("temperature", "density")

wt = "WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030"
@requires_ds(wt)
def test_wind_tunnel():
    ds = data_dir_load(wt)
    yield assert_equal, str(ds), "windtunnel_4lev_hdf5_plt_cnt_0030"
    for test in small_patch_amr(ds, _fields_2d):
        test_wind_tunnel.__name__ = test.description
        yield test

@requires_file(wt)
def test_FLASHDataset():
    assert isinstance(data_dir_load(wt), FLASHDataset)

@requires_file(sloshing)
def test_units_override():
    for test in units_override_check(sloshing):
        yield test

fid_1to3_b1 = "fiducial_1to3_b1/fiducial_1to3_b1_hdf5_part_0080"

fid_1to3_b1_fields = OrderedDict(
    [
        (("deposit", "all_density"), None),
        (("deposit", "all_count"), None),
        (("deposit", "all_cic"), None),
        (("deposit", "all_cic_velocity_x"), ("deposit", "all_cic")),
        (("deposit", "all_cic_velocity_y"), ("deposit", "all_cic")),
        (("deposit", "all_cic_velocity_z"), ("deposit", "all_cic")),
    ]
)


@requires_file(fid_1to3_b1)
def test_FLASHParticleDataset():
    assert isinstance(data_dir_load(fid_1to3_b1), FLASHParticleDataset)


dens_turb_mag = 'DensTurbMag/DensTurbMag_hdf5_plt_cnt_0015'
@requires_file(dens_turb_mag)
def test_FLASH25_dataset():
    ds = data_dir_load(dens_turb_mag)
    assert_equal(ds.parameters['time'], 751000000000.0)
    assert_equal(ds.domain_dimensions, np.array([8, 8, 8]))
    assert_equal(ds.domain_left_edge, 
                 ds.arr([-2e18, -2e18, -2e18], 'code_length'))

    assert_equal(ds.index.num_grids, 73)
    dd = ds.all_data()
    dd['density']


@requires_ds(fid_1to3_b1, big_data=True)
def test_fid_1to3_b1():
    ds = data_dir_load(fid_1to3_b1)
    for test in sph_answer(ds, 'fiducial_1to3_b1_hdf5_part_0080', 6684119, fid_1to3_b1_fields):
        test_fid_1to3_b1.__name__ = test.description
        yield test
