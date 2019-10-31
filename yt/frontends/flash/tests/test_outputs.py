"""
Title: test_flash.py
Purpose: FLASH frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from collections import OrderedDict

import numpy as np
import pytest

from yt.frontends.flash.api import FLASHDataset, \
    FLASHParticleDataset
from yt.testing import \
    assert_equal, \
    requires_file, \
    units_override_check
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

# Test data
sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
wt = "WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030"
fid_1to3_b1 = "fiducial_1to3_b1/fiducial_1to3_b1_hdf5_part_0080"
dens_turb_mag = 'DensTurbMag/DensTurbMag_hdf5_plt_cnt_0015'


#============================================
#                 TestFlash
#============================================
@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestFlash(fw.AnswerTest):
    #-----
    # test_sloshing
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(sloshing)
    def test_sloshing(self, a, d, w, f, ds_sloshing):
        # Run the small_patch_amr test suite
        self.hashes.update(self.small_patch_amr(ds_sloshing, f, w, a, d))

    #-----
    # test_wind_tunnel
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(wt)
    def test_wind_tunnel(self, a, d, w, f, ds_wt):
        # Run the small_patch_amr test suite
        self.hashes.update(self.small_patch_amr(ds_wt, f, w, a, d))

    #-----
    # test_FLASHDataset
    #-----
    @requires_file(wt)
    def test_FLASHDataset(self, ds_wt):
        assert isinstance(ds_wt, FLASHDataset)

    #-----
    # test_units_override
    #-----
    @requires_file(sloshing)
    def test_units_override(self, ds_sloshing):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        units_override_check(ds_sloshing, sloshing)

    #-----
    # test_FLASHParticleDataset
    #-----
    @requires_file(fid_1to3_b1)
    def test_FLASHParticleDataset(self, ds_fid_1to3_b1):
        assert isinstance(ds_fid_1to3_b1, FLASHParticleDataset)

    #-----
    # test_Flash25_dataset
    #-----
    @requires_file(dens_turb_mag)
    def test_FLASH25_dataset(self, ds_dens_turb_mag):
        assert_equal(ds_dens_turb_mag.parameters['time'], 751000000000.0)
        assert_equal(ds_dens_turb_mag.domain_dimensions, np.array([8, 8, 8]))
        assert_equal(ds_dens_turb_mag.domain_left_edge,
            ds_dens_turb_mag.arr([-2e18, -2e18, -2e18], 'code_length'))
        assert_equal(ds_dens_turb_mag.index.num_grids, 73)
        dd = ds_dens_turb_mag.all_data()
        dd['density']

    #-----
    # test_fid_1to3_b1
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(fid_1to3_b1)
    def test_fid_1to3_b1(self, f, w, d, a, ds_fid_1to3_b1):
        self.hashes.update(self.sph_answer(ds_fid_1to3_b1,
            'fiducial_1to3_b1_hdf5_part_0080', 6684119, f, w, d, a
            )
        )
