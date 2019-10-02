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
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('answer_file')
class TestGamer(fw.AnswerTest):
    #-----
    # test_jet
    #-----
    @pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
        reason="Skipping test_jet because --answer-big-data was not set."
    )
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(jet)
    def test_jet(self, ds_jet):
        ds = ds_jet
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        # Run the small_patch_amr test suite
        self.hashes = self.small_patch_amr(ds, _fields_jet, weights, axes, ds_objs)

    #-----
    # test_psiDM
    #-----
    @pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
        reason="Skipping test_jet because --answer-big-data was not set."
    )
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(psiDM)
    def test_psiDM(self, ds_psiDM):
        ds = ds_psiDM
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        # Run the small_patch_amr test suite
        self.hashes = self.small_patch_amr(ds, _fields_psiDM, weights, axes, ds_objs)

    #-----
    # test_plummer
    #-----
    @pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
        reason="Skipping test_jet because --answer-big-data was not set."
    )
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(plummer)
    def test_plummer(self):
        ds = utils.data_dir_load(plummer)
        assert_equal(str(ds), "plummer_000000")
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        # Run the small_patch_amr test suite
        self.hashes = self.small_patch_amr(ds, _fields_plummer, weights, axes, ds_objs)

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
