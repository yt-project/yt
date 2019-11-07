"""
Title: test_exodusii.py
Purpose: Exodus II frontend tests
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.testing import \
    assert_array_equal, \
    assert_equal, \
    requires_file 
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

# Test data
out = "ExodusII/out.e"
out_s002 = "ExodusII/out.e-s002"
gold = "ExodusII/gold.e"
big_data = "MOOSE_sample_data/mps_out.e"


#============================================
#               TestExodusII
#============================================
@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestExodusII(fw.AnswerTest):
    #-----
    # test_out
    #-----
    @requires_file(out)
    def test_out(self, ds_out):
        field_list = [('all', 'conv_indicator'),
                      ('all', 'conv_marker'),
                      ('all', 'convected'),
                      ('all', 'diffused'),
                      ('connect1', 'conv_indicator'),
                      ('connect1', 'conv_marker'),
                      ('connect1', 'convected'),
                      ('connect1', 'diffused'),
                      ('connect2', 'conv_indicator'),
                      ('connect2', 'conv_marker'),
                      ('connect2', 'convected'),
                      ('connect2', 'diffused')]
        assert_equal(ds_out.dimensionality, 3)
        assert_equal(ds_out.current_time, 0.0)
        assert_array_equal(ds_out.parameters['nod_names'], ['convected', 'diffused'])
        assert_equal(ds_out.parameters['num_meshes'], 2)
        assert_array_equal(ds_out.field_list, field_list)

    #-----
    # test_out002
    #-----
    @requires_file(out_s002)
    def test_out002(self, ds_out_s002):
        field_list = [('all', 'conv_indicator'),
                      ('all', 'conv_marker'),
                      ('all', 'convected'),
                      ('all', 'diffused'),
                      ('connect1', 'conv_indicator'),
                      ('connect1', 'conv_marker'),
                      ('connect1', 'convected'),
                      ('connect1', 'diffused'),
                      ('connect2', 'conv_indicator'),
                      ('connect2', 'conv_marker'),
                      ('connect2', 'convected'),
                      ('connect2', 'diffused')]
        assert_equal(ds_out_s002.dimensionality, 3)
        assert_equal(ds_out_s002.current_time, 2.0)
        assert_array_equal(ds_out_s002.field_list, field_list)

    #-----
    # test_gold
    #-----
    @requires_file(gold)
    def test_gold(self, ds_gold):
        field_list = [('all', 'forced'), ('connect1', 'forced')]
        assert_array_equal(ds_gold.field_list, field_list)

    #-----
    # test_displacement_fields
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(big_data)
    def test_displacement_fields(self, disp):
        self.hashes.update({'generic_array' : {}})
        ds = utils.data_dir_load(big_data, kwargs={'displacements':disp})
        for mesh in ds.index.meshes:
            def array_func(*args, **kwargs):
                return mesh.connectivity_coords
            ga_hd = self.generic_array_test(array_func, args=[12])
            self.hashes.update({'generic_array' :  {str(mesh) : ga_hd}})
