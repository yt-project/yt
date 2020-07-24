"""
Title: test_exodusii.py
Purpose: Exodus II frontend tests
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.testing import assert_array_equal, assert_equal
from yt.utilities.answer_testing.answer_tests import generic_array
from yt.utilities.answer_testing.utils import data_dir_load, requires_ds

# Test data
out = "ExodusII/out.e"
out_s002 = "ExodusII/out.e-s002"
gold = "ExodusII/gold.e"
big_data = "MOOSE_sample_data/mps_out.e"


@pytest.mark.answer_test
@pytest.mark.usefixtures("answer_file")
class TestExodusII:
    @pytest.mark.parametrize("ds", [out], indirect=True)
    def test_out(self, ds):
        field_list = [
            ("all", "conv_indicator"),
            ("all", "conv_marker"),
            ("all", "convected"),
            ("all", "diffused"),
            ("connect1", "conv_indicator"),
            ("connect1", "conv_marker"),
            ("connect1", "convected"),
            ("connect1", "diffused"),
            ("connect2", "conv_indicator"),
            ("connect2", "conv_marker"),
            ("connect2", "convected"),
            ("connect2", "diffused"),
        ]
        assert_equal(ds.dimensionality, 3)
        assert_equal(ds.current_time, 0.0)
        assert_array_equal(ds.parameters["nod_names"], ["convected", "diffused"])
        assert_equal(ds.parameters["num_meshes"], 2)
        assert_array_equal(ds.field_list, field_list)

    @pytest.mark.parametrize("ds", [out_s002], indirect=True)
    def test_out002(self, ds):
        field_list = [
            ("all", "conv_indicator"),
            ("all", "conv_marker"),
            ("all", "convected"),
            ("all", "diffused"),
            ("connect1", "conv_indicator"),
            ("connect1", "conv_marker"),
            ("connect1", "convected"),
            ("connect1", "diffused"),
            ("connect2", "conv_indicator"),
            ("connect2", "conv_marker"),
            ("connect2", "convected"),
            ("connect2", "diffused"),
        ]
        assert_equal(ds.dimensionality, 3)
        assert_equal(ds.current_time, 2.0)
        assert_array_equal(ds.field_list, field_list)

    @pytest.mark.parametrize("ds", [gold], indirect=True)
    def test_gold(self, ds):
        field_list = [("all", "forced"), ("connect1", "forced")]
        assert_array_equal(ds.field_list, field_list)

    @pytest.mark.usefixtures("hashing")
    @requires_ds(big_data)
    def test_displacement_fields(self, disp):
        self.hashes.update({"generic_array": {}})
        ds = data_dir_load(big_data, kwargs={"displacements": disp})
        for mesh in ds.index.meshes:

            def array_func(*args, **kwargs):
                return mesh.connectivity_coords

            ga = generic_array(array_func, args=[12])
            self.hashes.update({"generic_array": {str(mesh): ga}})
