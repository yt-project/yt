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
import sys
import tempfile
import shutil
from yt.testing import \
    fake_random_pf
from yt.mods import \
    SlicePlot, ProjectionPlot, OffAxisSlicePlot, OffAxisProjectionPlot


EXT_TO_TYPE = {
    '.ps': 'PostScript document text conforming DSC level 3.0',
    '.eps': 'PostScript document text conforming DSC level 3.0, type EPS',
    '.pdf': 'PDF document, version 1.4',
    '.png': 'PNG image data, 1070 x 1000, 8-bit/color RGBA, non-interlaced'
}


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def assert_fname(fname):
    """Function that checks file type using libmagic"""
    if fname is None:
        return

    try:
        import magic
    except ImportError:
        # OS X doesn't come with libmagic
        pass

    if 'magic' in sys.modules:
        ext = os.path.splitext(fname)[1]
        mds = magic.open(magic.MAGIC_NONE)
        mds.load()
        magic_text = mds.file(fname)
        mds.close()
        assert magic_text == EXT_TO_TYPE[ext]


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
            yield assert_fname, obj.save(fname)[0]

            obj = ProjectionPlot(test_pf, dim, 'Density')
            yield assert_fname, obj.save(fname)[0]

        obj = OffAxisSlicePlot(test_pf, normal, 'Density')
        yield assert_fname, obj.save(fname)[0]

        obj = OffAxisProjectionPlot(test_pf, normal, 'Density')
        yield assert_fname, obj.save(fname)[0]

    os.chdir(curdir)
    # clean up
    shutil.rmtree(tmpdir)
