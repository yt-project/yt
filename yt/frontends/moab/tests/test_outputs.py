"""
Title: test_moab.py
Purpose: Tests of semi-structured meshes in MoabHex8 format.
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from collections import OrderedDict

import numpy as np
import pytest

from yt.testing import \
    assert_equal, \
    assert_almost_equal, \
    requires_file, \
    units_override_check
from yt.frontends.moab.api import MoabHex8Dataset
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
c5 = "c5/c5.h5m"


# Globals
_fields = (("moab", "flux"),
          )


#============================================
#                  TestMoab
#============================================
@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestMoab(fw.AnswerTest):
    #-----
    # test_cantor_5
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(c5)
    def test_cantor_5(self, f, d, ds_c5):
        fv_hd = self.field_values_test(ds_c5, f, d)
        self.hashes.update({'field_values': fv_hd})

    @utils.requires_ds(c5)
    def test_cantor_5_non_param(self, ds_c5):
        ds = ds_c5
        np.random.seed(0x4d3d3d3)
        dd = ds.all_data()
        assert_almost_equal(ds.index.get_smallest_dx(), 0.00411522633744843, 10)
        assert_equal(dd["x"].shape[0], 63*63*63)
        assert_almost_equal(
            dd["cell_volume"].in_units("code_length**3").sum(dtype="float64").d,
            1.0, 10)
        for offset_1 in [1e-9, 1e-4, 0.1]:
            for offset_2 in [1e-9, 1e-4, 0.1]:
                DLE = ds.domain_left_edge
                DRE = ds.domain_right_edge
                ray = ds.ray(DLE + offset_1 * DLE.uq,
                             DRE - offset_2 * DRE.uq)
                assert_almost_equal(ray["dts"].sum(dtype="float64"), 1.0, 8)
        for i, p1 in enumerate(np.random.random((5, 3))):
            for j, p2 in enumerate(np.random.random((5, 3))):
                ray = ds.ray(p1, p2)
                assert_almost_equal(ray["dts"].sum(dtype="float64"), 1.0, 8)

    #-----
    # test_MoabHex8Dataset
    #-----
    @requires_file(c5)
    def test_MoabHex8Dataset(self, ds_c5):
        assert isinstance(ds_c5, MoabHex8Dataset)

    #-----
    # test_units_override
    #-----
    @requires_file(c5)
    def test_units_override(self, ds_c5):
        units_override_check(ds_c5, c5)
