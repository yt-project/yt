"""
ART frontend tests using SFG1 a=0.330

Author: Christopher Erick Moody <chrisemoody@gmail.com>
Affiliation: University of California Santa Cruz
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
from yt.frontends.art.api import ARTStaticOutput

_fields = ("Density", "particle_mass", ("all", "particle_position_x"))

sfg1 = "10MpcBox_csf512_a0.330.d"


@requires_pf(sfg1, big_data=True)
def test_sfg1():
    pf = data_dir_load(sfg1)
    yield assert_equal, str(pf), "10MpcBox_csf512_a0.330.d"
    dso = [None, ("sphere", ("max", (0.1, 'unitary')))]
    for field in _fields:
        for axis in [0, 1, 2]:
            for ds in dso:
                for weight_field in [None, "Density"]:
                    yield PixelizedProjectionValuesTest(
                        sfg1, axis, field, weight_field,
                        ds)
