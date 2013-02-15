"""
Orion frontend tests

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
from yt.frontends.orion.api import OrionStaticOutput

_fields = ("Temperature", "Density", "VelocityMagnitude", "DivV")

radadvect = "RadAdvect/plt00000"
@requires_pf(radadvect)
def test_radadvect():
    pf = data_dir_load(radadvect)
    yield assert_equal, str(pf), "plt00000"
    for test in small_patch_amr(radadvect, _fields):
        yield test

rt = "RadTube/plt00500"
@requires_pf(rt)
def test_radtube():
    pf = data_dir_load(rt)
    yield assert_equal, str(pf), "plt00500"
    for test in small_patch_amr(rt, _fields):
        yield test
