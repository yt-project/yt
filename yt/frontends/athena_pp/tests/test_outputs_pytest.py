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
from yt.testing import assert_allclose, requires_file, units_override_check
from yt.utilities.answer_testing import utils
from yt.utilities.answer_testing.answer_tests import field_values 
from yt.utilities.answer_testing.answer_tests import generic_array
from yt.utilities.answer_testing.answer_tests import grid_hierarchy
from yt.utilities.answer_testing.answer_tests import parentage_relationships
from yt.utilities.answer_testing.answer_tests import grid_values 
from yt.utilities.answer_testing.answer_tests import projection_values 

# Test data
disk = "KeplerianDisk/disk.out1.00000.athdf"
AM06 = "AM06/AM06.out1.00400.athdf"


uo_AM06 = {
    "length_unit": (1.0, "kpc"),
    "mass_unit": (1.0, "Msun"),
    "time_unit": (1.0, "Myr"),
}
AM06_kwargs = {"kwargs" : {"units_override" : uo_AM06}}
a_list = [0, 1, 2]
d_list = [None, ("sphere", ("max", (0.1, "unitary")))]
w_list = [None, "density"]
f_list = ["temperature", "density", "velocity_magnitude", "magnetic_field_x"]


@pytest.mark.answer_test
class TestAthenaPP:
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [disk], indirect=True)
    @pytest.mark.parametrize("f", ["density", "velocity_r"], indirect=True)
    def test_disk(self, field, ds):
        dd = ds.all_data()
        vol = (ds.domain_right_edge[0] ** 3 - ds.domain_left_edge[0] ** 3) / 3.0
        vol *= np.cos(ds.domain_left_edge[1]) - np.cos(ds.domain_right_edge[1])
        vol *= ds.domain_right_edge[2].v - ds.domain_left_edge[2].v
        assert_allclose(dd.quantities.total_quantity("cell_volume"), vol)

        def field_func(name):
            return dd[field]

        ga = generic_array(field_func, args=[field])
        self.hashes.update({"generic_array": ga})

    @pytest.mark.parametrize("ds", [[AM06, units_override=uo_AM06]], indirect=True)
    def test_AM06_override(self, ds):
        r"""Verify that overriding units causes derived unit values to be
        updated. See issue #1259.
        """
        assert float(ds.magnetic_unit.in_units("gauss")) == 9.01735778342523e-08

    @pytest.mark.parametrize("ds", [AM06], indirect=True)
    def test_units_override(self, ds):
        units_override_check(AM06)

    @pytest.mark.parametrize("ds", [AM06], indirect=True)
    def test_AthenaPPDataset(self, ds):
        assert isinstance(ds, AthenaPPDataset)

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [AM06], indirect=True)
    def test_gh_pr(self, ds):
        self.hashes.update({"grid_hierarchy": grid_hierarchy(ds)})
        self.hashes.update({"parentage_relationships": parentage_relationships(ds)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [AM06], indirect=True)
    @pytest.mark.parametrize("f", f_list, indirect=True)
    def test_gv(self, f, ds):
        self.hashes.update({"grid_values": grid_values(ds, f)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [AM06], indirect=True)
    @pytest.mark.parametrize("f", f_list, indirect=True)
    @pytest.mark.parametrize("d", d_list, indirect=True)
    def test_fv(self, d, f, ds):
        self.hashes.update({"field_values": field_values(ds, f, d)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [AM06], indirect=True)
    @pytest.mark.parametrize("f", f_list, indirect=True)
    @pytest.mark.parametrize("d", d_list, indirect=True)
    @pytest.mark.parametrize("a", a_list, indirect=True)
    @pytest.mark.parametrize("w", w_list, indirect=True)
    def test_pv(self, a, d, w, f, ds):
        self.hashes.update({"projection_values": projection_values(ds, a, f, w, d)})
