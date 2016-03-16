"""
Exodus II frontend tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import \
    assert_equal, \
    assert_array_equal, \
    requires_file
from yt.utilities.answer_testing.framework import \
    data_dir_load, \
    requires_ds, \
    GenericArrayTest


out = "ExodusII/out.e"


@requires_file(out)
def test_out():
    ds = data_dir_load(out)
    field_list = [('connect1', 'conv_indicator'),
                  ('connect1', 'conv_marker'),
                  ('connect1', 'convected'),
                  ('connect1', 'diffused'),
                  ('connect2', 'conv_indicator'),
                  ('connect2', 'conv_marker'),
                  ('connect2', 'convected'),
                  ('connect2', 'diffused')]
    yield assert_equal, str(ds), "out.e"
    yield assert_equal, ds.dimensionality, 3
    yield assert_equal, ds.current_time, 0.0
    yield assert_array_equal, ds.parameters['nod_names'], ['convected', 'diffused']
    yield assert_equal, ds.parameters['num_meshes'], 2
    yield assert_array_equal, ds.field_list, field_list 

out_s002 = "ExodusII/out.e-s002"


@requires_file(out_s002)
def test_out002():
    ds = data_dir_load(out_s002)
    field_list = [('connect1', 'conv_indicator'),
                  ('connect1', 'conv_marker'),
                  ('connect1', 'convected'),
                  ('connect1', 'diffused'),
                  ('connect2', 'conv_indicator'),
                  ('connect2', 'conv_marker'),
                  ('connect2', 'convected'),
                  ('connect2', 'diffused')]
    yield assert_equal, str(ds), "out.e-s002"
    yield assert_equal, ds.dimensionality, 3
    yield assert_equal, ds.current_time, 2.0
    yield assert_array_equal, ds.field_list, field_list 

gold = "ExodusII/gold.e"


@requires_file(gold)
def test_gold():
    ds = data_dir_load(gold)
    field_list = [('connect1', 'forced')]
    yield assert_equal, str(ds), "gold.e"
    yield assert_array_equal, ds.field_list, field_list 

big_data = "MOOSE_sample_data/mps_out.e"


@requires_ds(big_data)
def test_displacement_fields():
    displacement_dicts =[{'connect2': (5.0, [0.0, 0.0, 0.0])},
                         {'connect1': (1.0, [1.0, 2.0, 3.0]), 
                          'connect2': (0.0, [0.0, 0.0, 0.0])}]
    for disp in displacement_dicts:
        ds = data_dir_load(big_data, displacements=disp)
        for mesh in ds.index.meshes:
            def array_func():
                return mesh.connectivity_coords
            yield GenericArrayTest(ds, array_func, 12)
