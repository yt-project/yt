"""
title: test_chombo.py
Purpose: Chombo frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from yt.frontends.chombo.api import \
    ChomboDataset, \
    Orion2Dataset, \
    PlutoDataset
from yt.testing import \
    assert_equal, \
    requires_file, \
    units_override_check

import framework as fw
import utils

# Test data
gc = "GaussianCloud/data.0077.3d.hdf5"
tb = "TurbBoxLowRes/data.0005.3d.hdf5"
iso = "IsothermalSphere/data.0000.3d.hdf5"
zp = "ZeldovichPancake/plt32.2d.hdf5"
kho = "KelvinHelmholtz/data.0004.hdf5"


#============================================
#                 TestChombo
#============================================
class TestChombo(fw.AnswerTest):
    """
    Container for chombo frontend tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_gc
    #-----
    @utils.requires_ds(gc)
    def test_gc(self):
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
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("density", "velocity_magnitude", "magnetic_field_x")
        ds = utils.data_dir_load(gc)
        assert_equal(str(ds),  "data.0077.3d.hdf5")
        # Run small patch amr test suite
        hashes = self.small_patch_amr(ds, fields, weights, axes, ds_objs)
        # Save or compare hashes
        utils.handle_hashes(self.save_dir, 'gaussiancloud', hashes, self.answer_store)

    #-----
    # test_tb
    #-----
    @utils.requires_ds(tb)
    def test_tb(self):
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
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("density", "velocity_magnitude", "magnetic_field_x")
        ds = utils.data_dir_load(tb)
        assert_equal(str(ds), "data.0005.3d.hdf5")
        # Run small patch amr test suite
        hashes = self.small_patch_amr(ds, fields, weights, axes, ds_objs)
        # Save or compare hashes
        utils.handle_hashes(self.save_dir, 'turboxlowres', hashes, self.answer_store)

    #-----
    # test_iso
    #-----
    @utils.requires_ds(iso)
    def test_iso(self):
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
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("density", "velocity_magnitude", "magnetic_field_x")
        ds = utils.data_dir_load(iso)
        assert_equal(str(ds), "data.0000.3d.hdf5")
        # Run small patch amr test suite
        hashes = self.small_patch_amr(ds, fields, weights, axes, ds_objs)
        # Save or compare hashes
        utils.handle_hashes(self.save_dir, 'isothermalsphere', hashes, self.answer_store)

    #-----
    # test_zp
    #-----
    @utils.requires_ds(zp)
    def test_zp(self):
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
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("rhs", "phi")
        ds = utils.data_dir_load(zp)
        assert_equal(str(ds), "plt32.2d.hdf5")
        # Run small patch amr test suite
        hashes = self.small_patch_amr(ds, fields, weights, axes, ds_objs)
        # Save or compare hashes
        utils.handle_hashes(self.save_dir, 'zeldovichpancake', hashes, self.answer_store)

    #-----
    # test_kho
    #-----
    @utils.requires_ds(kho)
    def test_kho(self):
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
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("density", "velocity_magnitude", "magnetic_field_x")
        ds = utils.data_dir_load(kho)
        assert_equal(str(ds), "data.0004.hdf5")
        # Run small patch amr test suite
        hashes = self.small_patch_amr(ds, fields, weights, axes, ds_objs)
        # Save or compare hashes
        utils.handle_hashes(self.save_dir, 'kelvinhelmholtz', self.answer_store)

    #-----
    # test_ChomboDataset
    #-----
    @requires_file(zp)
    def test_ChomboDataset(self):
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
        assert isinstance(utils.data_dir_load(zp), ChomboDataset)

    #-----
    # test_Orion2Dataset
    #-----
    @requires_file(gc)
    def test_Orion2Dataset(self):
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
        assert isinstance(utils.data_dir_load(gc), Orion2Dataset)

    #-----
    # test_PlutoDataset
    #-----
    @requires_file(kho)
    def test_PlutoDataset(self):
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
        assert isinstance(utils.data_dir_load(kho), PlutoDataset)

    #-----
    # test_units_override_zp
    #-----
    @requires_file(zp)
    def test_units_override_zp(self):
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
        units_override_check(zp)

    #-----
    # test_units_override_gc
    #-----
    @requires_file(gc)
    def test_units_override_gc(self):
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
        units_override_check(gc)

    #-----
    # test_units_override_kho
    #-----
    @requires_file(kho)
    def test_units_override_kho(self):
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
        units_override_check(kho)
