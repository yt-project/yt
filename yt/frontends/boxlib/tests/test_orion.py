"""
Orion frontend tests



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
from yt.frontends.orion.api import OrionStaticOutput

_fields = ("Temperature", "Density", "VelocityMagnitude", "DivV")

radadvect = "RadAdvect/plt00000"
@requires_pf(radadvect)
def test_radadvect():
    pf = data_dir_load(radadvect)
    yield assert_equal, str(pf), "plt00000"
    for test in small_patch_amr(radadvect, _fields):
        test_radadvect.__name__ = test.description
        yield test

rt = "RadTube/plt00500"
@requires_pf(rt)
def test_radtube():
    pf = data_dir_load(rt)
    yield assert_equal, str(pf), "plt00500"
    for test in small_patch_amr(rt, _fields):
        test_radtube.__name__ = test.description
        yield test
