"""
Testsuite for PlotWindow class

Author: Nathan Goldbaum <goldbaum@ucolick.org>
Affiliation: UCSC Astronomy
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Nathan Goldbaum.  All Rights Reserved.

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
import os
import tempfile
import shutil
from yt.testing import \
    fake_random_pf, assert_equal
from yt.mods import \
    SlicePlot, ProjectionPlot, OffAxisSlicePlot, OffAxisProjectionPlot


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def assert_fname(fname):
    """Function that checks file type using libmagic"""
    if fname is None:
        return

    with open(fname, 'rb') as fimg:
        data = fimg.read()
    data = str(data)
    image_type = ''

    # see http://www.w3.org/TR/PNG/#5PNG-file-signature
    if data.startswith('\211PNG\r\n\032\n'):
        image_type = '.png'
    # see http://www.mathguide.de/info/tools/media-types/image/jpeg
    elif data.startswith('\377\330'):
        image_type = '.jpeg'
    elif data.startswith('%!PS-Adobe'):
        if 'EPSF' in data[:data.index('\n')]:
            image_type = '.eps'
        else:
            image_type = '.ps'
    elif data.startswith('%PDF'):
        image_type = '.pdf'

    return image_type == os.path.splitext(fname)[1]


def test_plotwindow():
    """Main test suite for PlotWindow."""
    # Perform I/O in safe place instead of yt main dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    normal = [1, 1, 1]

    test_pf = fake_random_pf(64)
    test_flnms = [None, 'test.png', 'test.eps',
                  'test.ps', 'test.pdf']
    for fname in test_flnms:
        for dim in [0, 1, 2]:
            obj = SlicePlot(test_pf, dim, 'Density')
            yield assert_equal, assert_fname(obj.save(fname)[0]), True

            obj = ProjectionPlot(test_pf, dim, 'Density')
            yield assert_equal, assert_fname(obj.save(fname)[0]), True

        obj = OffAxisSlicePlot(test_pf, normal, 'Density')
        yield assert_equal, assert_fname(obj.save(fname)[0]), True

        obj = OffAxisProjectionPlot(test_pf, normal, 'Density')
        yield assert_equal, assert_fname(obj.save(fname)[0]), True

    os.chdir(curdir)
    # clean up
    shutil.rmtree(tmpdir)
