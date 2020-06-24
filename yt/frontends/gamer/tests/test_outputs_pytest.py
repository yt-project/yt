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
    requires_file, \
    units_override_check
from yt.frontends.gamer.api import GAMERDataset
from yt.utilities.answer_testing.answer_tests import small_patch_amr
from yt.utilities.answer_testing import utils


# Test data
jet         = "InteractingJets/jet_000002"
psiDM       = "WaveDarkMatter/psiDM_000020"
plummer     = "Plummer/plummer_000000"
mhd_vortex   = "MHDOrszagTangVortex/Data_000018"


#============================================
#                 TestGamer
#============================================
@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestGamer:
    #-----
    # test_jet
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(jet)
    def test_jet(self, f, a, d, w, ds_jet):
        self.hashes.update(small_patch_amr(ds_jet, f, w, a, d))

    #-----
    # test_psiDM
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [psiDM], indirect=True)
    def test_psiDM(self, f, a, d, w, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    #-----
    # test_plummer
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [plummer], indirect=True)
    def test_plummer(self, f, a, d, w, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    #-----
    # test_mhd_vortex
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [mhd_vortex], indirect=True)
    def test_mhdvortex(self, f, a, d, w, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    #-----
    # test_GAMERDataset
    #-----
    @pytest.mark.parametrize('ds', [psiDM], indirect=True)
    def test_GAMERDataset(self, ds):
        assert isinstance(ds, GAMERDataset)

    #-----
    # test_units_override
    #-----
    @requires_file(jet)
    def test_units_override(self):
        units_override_check(jet)
