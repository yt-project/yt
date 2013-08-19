"""
Test for off_axis_projection and write_projection

Author: Cameron Hummels <chummels@gmail.com>
Affiliation: University of Arizona
Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Cameron Hummels.  All Rights Reserved.

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
import os.path
import tempfile
import shutil
from yt.testing import \
    fake_random_pf, assert_equal, expand_keywords
from yt.mods import \
    off_axis_projection, write_projection


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_write_projection():
    """Tests functionality of off_axis_projection and write_projection."""
    # Perform I/O in safe place instead of yt main dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    # args for off_axis_projection
    test_pf = fake_random_pf(64)
    c = [0.5, 0.5, 0.5]
    norm = [0.5, 0.5, 0.5]
    W = [0.5,0.5,1.0]
    N = 64
    field = "Density"
    oap_args = [test_pf, c, norm, W, N, field]

    # kwargs for off_axis_projection
    oap_kwargs = {}
    oap_kwargs['weight'] = (None, 'CellMassMsun')
    oap_kwargs['no_ghost'] = (True, False)
    oap_kwargs['interpolated'] = (True, False)
    oap_kwargs['north_vector'] = ((1,0,0), (0,0.5,1.0))
    oap_kwargs_list = expand_keywords(oap_kwargs)

    # args for write_projection
    fn = "test.png"

    # kwargs for write_projection
    wp_kwargs = {}
    wp_kwargs['colorbar'] = (True, False)
    wp_kwargs['colorbar_label'] = ('test')
    wp_kwargs['title'] = ('test')
    wp_kwargs['limits'] = (None, (1e3, 1e5))
    wp_kwargs['take_log'] = (True, False)
    wp_kwargs['figsize'] = ((8,6), [1,1])
    wp_kwargs['dpi'] = (100, 50)
    wp_kwargs['cmap_name'] = ('algae', 'jet')
    wp_kwargs_list = expand_keywords(wp_kwargs)

    # test all off_axis_projection kwargs and write_projection kwargs
    # make sure they are able to be projected, then remove and try next
    # iteration
    for oap_kwargs in oap_kwargs_list:
        image = off_axis_projection(*oap_args, **oap_kwargs)
        for wp_kwargs in wp_kwargs_list:
            write_projection(image, fn, **wp_kwargs)
            yield assert_equal, os.path.exists(fn), True
            os.remove(fn)

    os.chdir(curdir)
    # clean up
    shutil.rmtree(tmpdir)
