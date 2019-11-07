"""
Title: test_gamer.py
Purpose: GAMER frontend tests
Notes:
    Copyright (c) 2016, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.testing import \
    assert_equal, \
    requires_file, \
    units_override_check
from yt.frontends.gamer.api import GAMERDataset
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
jet         = "InteractingJets/jet_000002"
psiDM       = "WaveDarkMatter/psiDM_000020"
plummer     = "Plummer/plummer_000000"


# Globals
_fields_jet = ("temperature", "density", "velocity_magnitude")
_fields_psiDM = ("Dens", "Real", "Imag")
_fields_plummer = ( ("gamer","ParDens"), ("deposit","io_cic") )
jet_units   = {"length_unit":(1.0,"kpc"),
               "time_unit"  :(3.08567758096e+13,"s"),
               "mass_unit"  :(1.4690033e+36,"g")}


#============================================
#                 TestGamer
#============================================
@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestGamer(fw.AnswerTest):
    #-----
    # test_jet
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(jet)
    def test_jet(self, f, a, d, w, ds_jet):
        # Run the small_patch_amr test suite
        self.hashes.update(self.small_patch_amr(ds_jet, f, w, a, d))

    #-----
    # test_psiDM
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(psiDM)
    def test_psiDM(self, f, a, d, w, ds_psiDM):
        # Run the small_patch_amr test suite
        self.hashes.update(self.small_patch_amr(ds_psiDM, f, w, a, d))

    #-----
    # test_plummer
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(plummer)
    def test_plummer(self, f, a, d, w, ds_plummer):
        # Run the small_patch_amr test suite
        self.hashes.update(self.small_patch_amr(ds_plummer, f, w, a, d))

    #-----
    # test_mhd_vortex
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(mhd_vortex)
    def test_mhdvortex(self, f, a, d, w, ds_mhd_vortex):
        self.hashes.update(self.small_patch_amr(ds_mhd_vortex, f, w, a, d))

    #-----
    # test_GAMERDataset
    #-----
    @requires_file(psiDM)
    def test_GAMERDataset(self, ds_psiDM):
        assert isinstance(ds_psiDM, GAMERDataset)

    #-----
    # test_units_override
    #-----
    @requires_file(jet)
    def test_units_override(self, ds_jet):
        units_override_check(ds_jet, jet)
