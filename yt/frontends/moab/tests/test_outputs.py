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


# Answer file
answer_file = 'moab_answers.yaml'


#============================================
#                  TestMoab
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestMoab(fw.AnswerTest):
    """
    Container for moab frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_cantor_5
    #-----
    @utils.requires_ds(c5)
    def test_cantor_5(self, ds_c5):
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
        ds = ds_c5
        hashes = OrderedDict()
        hashes['field_values'] = OrderedDict()
        np.random.seed(0x4d3d3d3)
        dso = [ None, ("sphere", ("c", (0.1, 'unitary'))),
                      ("sphere", ("c", (0.2, 'unitary')))]
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
        for field in _fields:
            hashes['field_values'][field] = OrderedDict()
            for dobj_name in dso:
                fv_hd = utils.generate_hash(
                    self.field_values_test(ds, field, dobj_name)
                )
                hashes['field_values'][field][dobj_name] = fv_hd
        hashes = {'cantor_5' : hashes}
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

    #-----
    # test_MoabHex8Dataset
    #-----
    @requires_file(c5)
    def test_MoabHex8Dataset(self, ds_c5):
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
        assert isinstance(ds_c5, MoabHex8Dataset)

    #-----
    # test_units_override
    #-----
    @requires_file(c5)
    def test_units_override(self, ds_c5):
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
        units_override_check(ds_c5, c5)
