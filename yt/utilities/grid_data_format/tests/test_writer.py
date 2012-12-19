"""
Testsuite for writing yt data to GDF

Author: Kacper Kowalik <xarthisius.kk@gmail.com>
Affiliation: Torun Center for Astronomy, NCU
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Kacper Kowalik.  All Rights Reserved.

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
import tempfile
import shutil
import os
import h5py as h5
from yt.testing import \
    fake_random_pf, assert_equal
from yt.utilities.grid_data_format.writer import \
    write_to_gdf
from yt.frontends.gdf.data_structures import \
    GDFStaticOutput
from yt.mods import \
    load

TEST_AUTHOR = "yt test runner"
TEST_COMMENT = "Testing write_to_gdf"


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_write_gdf():
    """Main test suite for write_gdf"""
    tmpdir = tempfile.mkdtemp()
    tmpfile = os.path.join(tmpdir, 'test_gdf.h5')

    test_pf = fake_random_pf(64)
    write_to_gdf(test_pf, tmpfile, data_author=TEST_AUTHOR,
                 data_comment=TEST_COMMENT)
    del test_pf

    assert isinstance(load(tmpfile), GDFStaticOutput)

    h5f = h5.File(tmpfile, 'r')
    gdf = h5f['gridded_data_format'].attrs
    assert_equal(gdf['data_author'], TEST_AUTHOR)
    assert_equal(gdf['data_comment'], TEST_COMMENT)
    h5f.close()

    shutil.rmtree(tmpdir)
