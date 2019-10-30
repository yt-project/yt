"""
Title: test_fits.py
Purpose: FITS frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.testing import \
    assert_equal, \
    requires_file, \
    units_override_check
from yt.frontends.fits.data_structures import FITSDataset, \
    SpectralCubeFITSDataset, \
    SkyDataFITSDataset, \
    EventsFITSDataset
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
grs = "radio_fits/grs-50-cube.fits"
vf = "UnigridData/velocity_field_20.fits"
acis = "xray_fits/acisf05356N003_evt2.fits"
A2052 = "xray_fits/A2052_merged_0.3-2_match-core_tmap_bgecorr.fits"


# Globals
_fields_grs = ("temperature",)
_fields_vels = ("velocity_x","velocity_y","velocity_z")
_fields_acis = ("counts_0.1-2.0", "counts_2.0-5.0")
_fields_A2052 = ("flux",)


#============================================
#                 TestFits
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason=("--with-answer-testing not set."))
@pytest.mark.usefixtures('answer_file')
class TestFits(fw.AnswerTest):
    #-----
    # test_grs
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(grs)
    def test_grs(self, f, a, d, w, ds_grs):
        # Run the small_patch_amr test suite
        self.hashes.update(self.small_patch_amr(ds_grs, f, w, a, d))

    #-----
    # test_velocity_field
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(vf)
    def test_velocity_field(self, f, a, d, w, ds_vf):
        # Run the small_patch_amr test suite
        self.hashes.update(self.small_patch_amr(ds_vf, f, w, a, d))

    #-----
    # test_acts
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(acis)
    def test_acis(self, f, a, d, w, ds_acis):
        # Run the small_patch_amr test suite
        self.hashes.update(self.small_patch_amr(ds_acis, f, w, a, d))

    #-----
    # test_A2052
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(A2052)
    def test_A2052(self, f, a, d, w, ds_A2052):
        # Run the small_patch_amr test suite
        self.hashes.update(self.small_patch_amr(ds_A2052, f, w, a, d))

    #-----
    # test_units_override
    #-----
    @requires_file(vf)
    def test_units_override(self, ds_vf):
        units_override_check(ds_vf, vf)

    #-----
    # test_FITSDataset
    #-----
    @requires_file(vf)
    def test_FITSDataset(self, ds_vf):
        assert isinstance(ds_vf, FITSDataset)

    #-----
    # test_SpectralCubeFITSDataset
    #-----
    @requires_file(grs)
    def test_SpectralCubeFITSDataset(self, ds_grs):
        assert isinstance(ds_grs, SpectralCubeFITSDataset)

    #-----
    # test_EventsFITSDataset
    #-----
    @requires_file(acis)
    def test_EventsFITSDataset(self, ds_acis):
        assert isinstance(ds_acis, EventsFITSDataset)

    #-----
    # test_SkyDataFITSDataset
    #-----
    @requires_file(A2052)
    def test_SkyDataFITSDataset(self, ds_A2052):
        assert isinstance(ds_A2052, SkyDataFITSDataset)
