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
    data_dir_load

out = "ExodusII/out.e"


@requires_file(out)
def test_out():
    ds = data_dir_load(out)
    yield assert_equal, str(ds), "out.e"
    yield assert_equal, ds.dimensionality, 3
    yield assert_equal, ds.unique_identifier, 5081193338833632556
    yield assert_equal, ds.current_time, 0.0
    yield assert_array_equal, ds.parameters['nod_names'], ['convected', 'diffused']
    yield assert_equal, ds.parameters['num_meshes'], 2

out_s002 = "ExodusII/out.e-s002"


@requires_file(out_s002)
def test_out002():
    ds = data_dir_load(out_s002)
    yield assert_equal, str(ds), "out.e-s002"
    yield assert_equal, ds.dimensionality, 3
    yield assert_equal, ds.current_time, 2.0

gold = "ExodusII/gold.e"


@requires_file(gold)
def test_gold():
    ds = data_dir_load(gold)
    yield assert_equal, str(ds), "gold.e"
