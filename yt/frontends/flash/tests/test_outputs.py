"""
FLASH frontend tests



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
from yt.frontends.flash.api import FLASHStaticOutput

_fields = ("Temperature", "Density", "VelocityMagnitude", "DivV")

sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
@requires_pf(sloshing, big_data=True)
def test_sloshing():
    pf = data_dir_load(sloshing)
    yield assert_equal, str(pf), "sloshing_low_res_hdf5_plt_cnt_0300"
    for test in small_patch_amr(sloshing, _fields):
        test_sloshing.__name__ = test.description
        yield test

_fields_2d = ("Temperature", "Density")

wt = "WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030"
@requires_pf(wt)
def test_wind_tunnel():
    pf = data_dir_load(wt)
    yield assert_equal, str(pf), "windtunnel_4lev_hdf5_plt_cnt_0030"
    for test in small_patch_amr(wt, _fields_2d):
        test_wind_tunnel.__name__ = test.description
        yield test
