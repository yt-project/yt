"""
Title: test_gdf.py
Purpose: GDF frontend tests
Notes:
    Copyright (c) 2016, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.frontends.gdf.api import GDFDataset
from yt.testing import \
    requires_file, \
    units_override_check
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

# Test data
sedov = "sedov/sedov_tst_0004.h5"


#============================================
#                  TestGDF
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('answer_test')
class TestGDF(fw.AnswerTest):
    #-----
    # test_sedov_tunnel
    #-----
    @pytest.mark.usefixtures('hashing')
    @requires_file(sedov)
    def test_sedov_tunnel(self, ds_sedov):
        # Set up arrays for testing
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("density", "velocity_x")
        # Run the small_patch_amr test suite
        self.hashes = self.small_patch_amr(ds_sedov, fields, weights, axes, ds_objs)

    #-----
    # test_GDFDataset
    #-----
    @requires_file(sedov)
    def test_GDFDataset(self, ds_sedov):
        assert isinstance(ds_sedov, GDFDataset)

    #-----
    # test_units_override
    #-----
    @requires_file(sedov)
    def test_units_override(self, ds_sedov):
        units_override_check(ds_sedov, sedov)
