"""
FLASH frontend tests

Author: John ZuHone <jzuhone@gmail.com>
Affiliation: NASA/Goddard Space Flight Center
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 John ZuHone.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.testing import *
from yt.utilities.answer_testing.framework import \
    requires_pf, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load
from yt.frontends.flash.api import FLASHStaticOutput

_fields = ("Temperature", "Density", "VelocityMagnitude", "DivV")

sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
@requires_pf(sloshing)
def test_sloshing():
    pf = data_dir_load(sloshing)
    yield assert_equal, str(pf), "sloshing_low_res_hdf5_plt_cnt_0300"
    for test in small_patch_amr(sloshing, _fields):
        yield test

_fields_2d = ("Temperature", "Density")

wt = "WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030"
@requires_pf(wt)
def test_wind_tunnel():
    pf = data_dir_load(wt)
    yield assert_equal, str(pf), "windtunnel_4lev_hdf5_plt_cnt_0030"
    for test in small_patch_amr(wt, _fields_2d):
        yield test

gcm = "GalaxyClusterMerger/fiducial_1to10_b0.273d_hdf5_plt_cnt_0245.gz"
@requires_pf(gcm, big_data=True)
def test_galaxy_cluster_merger():
    pf = data_dir_load(gcm)
    for test in big_patch_amr(gcm, _fields):
        yield test

