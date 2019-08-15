"""
Title: test_fits.py
Purpose: FITS frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from yt.testing import \
    assert_equal, \
    requires_file, \
    units_override_check
from yt.frontends.fits.data_structures import FITSDataset, \
    SpectralCubeFITSDataset, \
    SkyDataFITSDataset, \
    EventsFITSDataset

import framework as fw
import utils


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
class TestFits(fw.AnswerTest):
    """
    Container for fits frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_grs
    #-----
    @utils.requires_ds(grs)
    def test_grs(self):
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
        ds = utils.data_dir_load(grs, cls=SpectralCubeFITSDataset, kwargs={"nan_mask":0.0})
        assert_equal(str(ds), "grs-50-cube.fits")
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "ones"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _fields_grs, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'test-fits-grs', hashes, self.answer_store)

    #-----
    # test_velocity_field
    #-----
    @utils.requires_ds(vf)
    def test_velocity_field(self):
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
        ds = utils.data_dir_load(vf, cls=FITSDataset)
        assert_equal(str(ds), "velocity_field_20.fits")
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "ones"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _fields_vels, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'test-fits-velocity', hashes, self.answer_store)

    #-----
    # test_acts
    #-----
    @utils.requires_ds(acis)
    def test_acis(self):
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
        from yt.frontends.fits.misc import setup_counts_fields
        ds = utils.data_dir_load(acis, cls=EventsFITSDataset)
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
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'test-fits-acis', hashes, self.answer_store)

    #-----
    # test_A2052
    #-----
    @utils.requires_ds(A2052)
    def test_A2052(self):
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
        ds = utils.data_dir_load(A2052, cls=SkyDataFITSDataset)
        assert_equal(str(ds), "A2052_merged_0.3-2_match-core_tmap_bgecorr.fits")
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "ones"]
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, _fields_A2052, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'test-fits-acis', hashes, self.answer_store)

    #-----
    # test_units_override
    #-----
    @requires_file(vf)
    def test_units_override(self):
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
        units_override_check(vf)

    #-----
    # test_FITSDataset
    #-----
    @requires_file(vf)
    def test_FITSDataset(self):
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
        assert isinstance(utils.data_dir_load(vf), FITSDataset)

    #-----
    # test_SpectralCubeFITSDataset
    #-----
    @requires_file(grs)
    def test_SpectralCubeFITSDataset(self):
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
        assert isinstance(utils.data_dir_load(grs), SpectralCubeFITSDataset)

    #-----
    # test_EventsFITSDataset
    #-----
    @requires_file(acis)
    def test_EventsFITSDataset(self):
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
        assert isinstance(utils.data_dir_load(acis), EventsFITSDataset)

    #-----
    # test_SkyDataFITSDataset
    #-----
    @requires_file(A2052)
    def test_SkyDataFITSDataset(self):
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
        assert isinstance(utils.data_dir_load(A2052), SkyDataFITSDataset)
