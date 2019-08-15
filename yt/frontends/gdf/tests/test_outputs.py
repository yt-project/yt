"""
Title: test_gdf.py
Purpose: GDF frontend tests
Notes:
    Copyright (c) 2016, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from yt.testing import units_override_check
from yt.frontends.gdf.api import GDFDataset

import framework as fw
import utils

# Test data
sedov = "sedov/sedov_tst_0004.h5"


#============================================
#                  TestGDF
#============================================
class TestGDF(fw.AnswerTest):
    """
    Conainer for gdf frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_sedov_tunnel
    #-----
    def test_sedov_tunnel(self):
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
        ds = utils.data_dir_load(sedov)
        assert str(ds) == "sedov_tst_0004"
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("density", "velocity_x")
        # Run the small_patch_amr test suite
        hashes = self.small_patch_amr(ds, fields, weights, axes, ds_objs)
        # Save or compare answer
        utils.handle_hashes(self.save_dir, 'sedov', hashes, self.answer_store)

    #-----
    # test_GDFDataset
    #-----
    def test_GDFDataset(self):
        assert isinstance(utils.data_dir_load(sedov), GDFDataset)

    #-----
    # test_units_override
    #-----
    def test_units_override(self):
        units_override_check(sedov)
