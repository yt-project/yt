"""
title: test_athena_pp.py
Purpose: Athena++ frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import numpy as np
import pytest

from yt.convenience import load
from yt.frontends.athena_pp.api import AthenaPPDataset
from yt.testing import \
    assert_allclose, \
    requires_file, \
    units_override_check
from yt.utilities.answer_testing.answer_tests import generic_array, \
    small_patch_amr
from yt.utilities.answer_testing import utils


# Test data
disk = "KeplerianDisk/disk.out1.00000.athdf"
AM06 = "AM06/AM06.out1.00400.athdf"


@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestAthenaPP:
    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [disk], indirect=True)
    def test_disk(self, field, ds):
        dd = ds.all_data()
        vol = (ds.domain_right_edge[0]**3-ds.domain_left_edge[0]**3)/3.0
        vol *= np.cos(ds.domain_left_edge[1])-np.cos(ds.domain_right_edge[1])
        vol *= ds.domain_right_edge[2].v-ds.domain_left_edge[2].v
        assert_allclose(dd.quantities.total_quantity("cell_volume"), vol)
        def field_func(name):
            return dd[field]
        ga = generic_array(field_func, args=[field])
        self.hashes.update({'generic_array' : ga})

    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [AM06], indirect=True)
    def test_AM06(self, a, d, w, f, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @requires_file(AM06)
    def test_AM06_override(self):
        r"""Verify that overriding units causes derived unit values to be
        updated. See issue #1259.
        """
        uo_AM06 = {
            'length_unit': (1.0, 'kpc'),
            'mass_unit': (1.0, 'Msun'),
            'time_unit': (1.0, 'Myr'),
        }
        ds = load(AM06, units_override=uo_AM06)
        assert float(ds.magnetic_unit.in_units('gauss')) == 9.01735778342523e-08

    @requires_file(AM06)
    def test_units_override(self):
        units_override_check(AM06)

    @requires_file(AM06)
    def test_AthenaPPDataset(self):
        assert isinstance(utils.data_dir_load(AM06), AthenaPPDataset)
