"""
ARTIO frontend tests 

Author: Samuel Leitner <sam.leitner@gmail.com>
Affiliation: University of Maryland College Park
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
    data_dir_load, \
    PixelizedProjectionValuesTest, \
    FieldValuesTest
from yt.frontends.artio.api import ARTIOStaticOutput

_fields = ("Temperature", "Density", "VelocityMagnitude",
           ("deposit", "all_density"), ("deposit", "all_count")) 

sizmbhloz = "sizmbhloz-clref04SNth-rs9_a0.9011/sizmbhloz-clref04SNth-rs9_a0.9011.art"
@requires_pf(sizmbhloz)
def test_sizmbhloz():
    pf = data_dir_load(sizmbhloz)
    yield assert_equal, str(pf), "sizmbhloz-clref04SNth-rs9_a0.9011.art"
    dso = [ None, ("sphere", ("max", (0.1, 'unitary')))]
    for field in _fields:
        for axis in [0, 1, 2]:
            for ds in dso:
                for weight_field in [None, "Density"]:
                    yield PixelizedProjectionValuesTest(
                        sizmbhloz, axis, field, weight_field,
                        ds)
                yield FieldValuesTest(
                        sizmbhloz, field, ds)
