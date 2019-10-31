"""
title: test_chombo.py
Purpose: Chombo frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.testing import \
    requires_file, \
    units_override_check
from yt.frontends.chombo.api import \
    ChomboDataset, \
    Orion2Dataset, \
    PlutoDataset
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

# Test data
gc = "GaussianCloud/data.0077.3d.hdf5"
tb = "TurbBoxLowRes/data.0005.3d.hdf5"
iso = "IsothermalSphere/data.0000.3d.hdf5"
zp = "ZeldovichPancake/plt32.2d.hdf5"
kho = "KelvinHelmholtz/data.0004.hdf5"


#============================================
#                 TestChombo
#============================================
@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestChombo(fw.AnswerTest):
    #-----
    # test_gc
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(gc)
    def test_gc(self, a, d, w, f, ds_gc):
        # Run small patch amr test suite
        self.hashes.update(self.small_patch_amr(ds_gc, f, w, a, d))

    #-----
    # test_tb
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(tb)
    def test_tb(self, a, d, w, f, ds_tb):
        # Run small patch amr test suite
        self.hashes.update(self.small_patch_amr(ds_tb, f, w, a, d))

    #-----
    # test_iso
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(iso)
    def test_iso(self, a, d, w, f, ds_iso):
        # Run small patch amr test suite
        self.hashes.update(self.small_patch_amr(ds_iso, f, w, a, d))

    #-----
    # test_zp
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(zp)
    def test_zp(self, a, d, w, f, ds_zp):
        # Run small patch amr test suite
        self.hashes.update(self.small_patch_amr(ds_zp, f, w, a, d))

    #-----
    # test_kho
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(kho)
    def test_kho(self, a, d, w, f, ds_kho):
        # Run small patch amr test suite
        self.hashes.update(self.small_patch_amr(ds_kho, f, w, a, d))

    #-----
    # test_ChomboDataset
    #-----
    @requires_file(zp)
    def test_ChomboDataset(self, ds_zp):
        assert isinstance(ds_zp, ChomboDataset)

    #-----
    # test_Orion2Dataset
    #-----
    @requires_file(gc)
    def test_Orion2Dataset(self, ds_gc):
        assert isinstance(ds_gc, Orion2Dataset)

    #-----
    # test_PlutoDataset
    #-----
    @requires_file(kho)
    def test_PlutoDataset(self, ds_kho):
        assert isinstance(ds_kho, PlutoDataset)

    #-----
    # test_units_override_zp
    #-----
    @requires_file(zp)
    def test_units_override_zp(self, ds_zp):
        units_override_check(ds_zp, zp)

    #-----
    # test_units_override_gc
    #-----
    @requires_file(gc)
    def test_units_override_gc(self, ds_gc):
        units_override_check(ds_gc, gc)

    #-----
    # test_units_override_kho
    #-----
    @requires_file(kho)
    def test_units_override_kho(self, ds_kho):
        units_override_check(ds_kho, kho)
