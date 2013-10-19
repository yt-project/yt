"""
Chombo frontend tests



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
    requires_pf, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load
from yt.frontends.chombo.api import ChomboStaticOutput

_fields = ("Density", "VelocityMagnitude", "DivV", "X-magnfield")

gc = "GaussianCloud/data.0077.3d.hdf5"
@requires_pf(gc)
def test_gc():
    pf = data_dir_load(gc)
    yield assert_equal, str(pf), "data.0077.3d.hdf5"
    for test in small_patch_amr(gc, _fields):
        test_gc.__name__ = test.description
        yield test

tb = "TurbBoxLowRes/data.0005.3d.hdf5"
@requires_pf(tb)
def test_tb():
    pf = data_dir_load(tb)
    yield assert_equal, str(pf), "data.0005.3d.hdf5"
    for test in small_patch_amr(tb, _fields):
        test_tb.__name__ = test.description
        yield test
