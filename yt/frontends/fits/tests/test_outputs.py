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


# Answer file
answer_file = 'fits_answers.yaml'


#============================================
#                 TestFits
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason=("--with-answer-testing not set."))
class TestFits(fw.AnswerTest):
    #-----
    # test_grs
    #-----
    @utils.requires_ds(grs)
    def test_grs(self, ds_grs):
        ds = ds_grs
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "ones"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _fields_grs, weights, axes, ds_objs)
        hashes = {'grs' : hashes}
        # Save or compare answer
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

    #-----
    # test_velocity_field
    #-----
    @utils.requires_ds(vf)
    def test_velocity_field(self, ds_vf):
        ds = ds_vf
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "ones"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _fields_vels, weights, axes, ds_objs)
        hashes = {'velocity_field' : hashes}
        # Save or compare answer
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

    #-----
    # test_acts
    #-----
    @utils.requires_ds(acis)
    def test_acis(self, ds_acis):
        from yt.frontends.fits.misc import setup_counts_fields
        ds = ds_acis
        ebounds = [(0.1, 2.0), (2.0, 5.0)]
        setup_counts_fields(ds, ebounds)
        assert_equal(str(ds), "acisf05356N003_evt2.fits")
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "ones"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _fields_acis, weights, axes, ds_objs)
        hashes = {'acis' : hashes}
        # Save or compare answer
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

    #-----
    # test_A2052
    #-----
    @utils.requires_ds(A2052)
    def test_A2052(self, ds_A2052):
        ds = ds_A2052
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "ones"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _fields_A2052, weights, axes, ds_objs)
        hashes = {'A2052' : hashes}
        # Save or compare answer
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

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
