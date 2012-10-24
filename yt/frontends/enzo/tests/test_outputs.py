"""
Enzo frontend tests using moving7

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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
from yt.frontends.enzo.api import EnzoStaticOutput

_fields = ("Temperature", "Density", "VelocityMagnitude", "DivV",
           "particle_density")

m7 = "DD0010/moving7_0010"
@requires_pf(m7)
def test_moving7():
    pf = data_dir_load(m7)
    yield assert_equal, str(pf), "moving7_0010"
    for test in small_patch_amr(m7, _fields):
        yield test

g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"
@requires_pf(g30, big_data=True)
def test_galaxy0030():
    pf = data_dir_load(g30)
    yield assert_equal, str(pf), "galaxy0030"
    for test in big_patch_amr(g30, _fields):
        yield test
