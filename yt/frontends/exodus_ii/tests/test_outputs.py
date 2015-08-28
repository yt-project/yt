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

from yt.testing import *
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load
from yt.frontends.exodus_ii.api import ExodusIIDataset
from yt.config import ytcfg
from yt.convenience import load

_fields = ("convected", "diffused")

out = "ExodusII/out.e"
@requires_ds(out)
def test_out():
    ds = data_dir_load(out)
    yield assert_equal, str(ds), "out.e"
    yield assert_equal, ds.dimensionality, 3
    yield assert_equal, ds.unique_identifier, 'Wed Apr 15 07:52:29 2015'
    yield assert_equal, ds.current_time, 0.0
    
out_s002 = "ExodusII/out.e-s002"
@requires_ds(out_s002)
def test_out002():
    ds = data_dir_load(out_s002)
    yield assert_equal, str(ds), "out.e-s002"
    yield assert_equal, ds.dimensionality, 3
    yield assert_equal, ds.current_time, 0.0
    


