"""
Chombo frontend tests

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Matthew Turk.  All Rights Reserved.

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
from yt.frontends.chombo.api import ChomboStaticOutput

_fields = ("Density", "VelocityMagnitude", "DivV", "X-magnfield")

gc = "GaussianCloud/data.0077.3d.hdf5"
@requires_pf(gc)
def test_gc():
    pf = data_dir_load(gc)
    yield assert_equal, str(pf), "data.0077.3d.hdf5"
    for test in small_patch_amr(gc, _fields):
        yield test

tb = "TurbBoxLowRes/data.0005.3d.hdf5"
@requires_pf(tb)
def test_tb():
    pf = data_dir_load(tb)
    yield assert_equal, str(pf), "data.0005.3d.hdf5"
    for test in small_patch_amr(tb, _fields):
        yield test
