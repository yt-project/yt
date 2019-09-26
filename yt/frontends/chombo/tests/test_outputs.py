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
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('answer_file')
class TestChombo(fw.AnswerTest):
    #-----
    # test_gc
    #-----
    @utils.requires_ds(gc)
    def test_gc(self, ds_gc):
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("density", "velocity_magnitude", "magnetic_field_x")
        # Run small patch amr test suite
        hashes = self.small_patch_amr(ds_gc, fields, weights, axes, ds_objs)
        hashes = {'gc' : hashes}
        # Save or compare hashes
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)

    #-----
    # test_tb
    #-----
    @utils.requires_ds(tb)
    def test_tb(self, ds_tb):
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("density", "velocity_magnitude", "magnetic_field_x")
        # Run small patch amr test suite
        hashes = self.small_patch_amr(ds_tb, fields, weights, axes, ds_objs)
        # Save or compare hashes
        hashes = {'tb' : hashes}
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)

    #-----
    # test_iso
    #-----
    @utils.requires_ds(iso)
    def test_iso(self, ds_iso):
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("density", "velocity_magnitude", "magnetic_field_x")
        # Run small patch amr test suite
        hashes = self.small_patch_amr(ds_iso, fields, weights, axes, ds_objs)
        hashes = {'iso' : hashes}
        # Save or compare hashes
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)

    #-----
    # test_zp
    #-----
    @utils.requires_ds(zp)
    def test_zp(self, ds_zp):
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "rhs"]
        fields = ("rhs", "phi")
        # Run small patch amr test suite
        hashes = self.small_patch_amr(ds_zp, fields, weights, axes, ds_objs)
        hashes = {'zp' : hashes}
        # Save or compare hashes
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)

    #-----
    # test_kho
    #-----
    @utils.requires_ds(kho)
    def test_kho(self, ds_kho):
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("density", "velocity_magnitude", "magnetic_field_x")
        # Run small patch amr test suite
        hashes = self.small_patch_amr(ds_kho, fields, weights, axes, ds_objs)
        hashes = {'kho' : hashes}
        # Save or compare hashes
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)

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
