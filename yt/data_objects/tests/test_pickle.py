"""
Testsuite for pickling yt objects.

Author: Elizabeth Tasker <tasker@astro1.sci.hokudai.ac.jp>
Affiliation: Hokkaido University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Elizabeth Tasker. All Rights Reserved.

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
import cPickle
import os
import tempfile
from yt.testing \
    import fake_random_pf, assert_equal


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_save_load_pickle():
    """Main test for loading pickled objects"""
    return # Until boolean regions are implemented we can't test this
    test_pf = fake_random_pf(64)

    # create extracted region from boolean (fairly complex object)
    center = (test_pf.domain_left_edge + test_pf.domain_right_edge) / 2
    sp_outer = test_pf.h.sphere(center, test_pf.domain_width[0])
    sp_inner = test_pf.h.sphere(center, test_pf.domain_width[0] / 10.0)
    sp_boolean = test_pf.h.boolean([sp_outer, "NOT", sp_inner])

    minv, maxv = sp_boolean.quantities["Extrema"]("Density")[0]
    contour_threshold = min(minv * 10.0, 0.9 * maxv)

    contours = sp_boolean.extract_connected_sets(
        "Density", 1, contour_threshold, maxv + 1, log_space=True, cache=True)

    # save object
    cpklfile = tempfile.NamedTemporaryFile(delete=False)
    cPickle.dump(contours[1][0], cpklfile)
    cpklfile.close()

    # load object
    test_load = cPickle.load(open(cpklfile.name, "rb"))

    assert_equal.description = \
        "%s: File was pickle-loaded succesfully" % __name__
    yield assert_equal, test_load is not None, True
    assert_equal.description = \
        "%s: Length of pickle-loaded connected set object" % __name__
    yield assert_equal, len(contours[1][0]), len(test_load)

    os.remove(cpklfile.name)
