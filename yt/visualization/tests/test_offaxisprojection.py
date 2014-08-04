"""
Test for off_axis_projection and write_projection



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import os
import os.path
import tempfile
import shutil
from yt.testing import \
    fake_random_ds, assert_equal, expand_keywords
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
    test_ds = fake_random_ds(64)
    c = [0.5, 0.5, 0.5]
    norm = [0.5, 0.5, 0.5]
    W = [0.5,0.5,1.0]
    N = 64
    field = "density"
    oap_args = [test_ds, c, norm, W, N, field]

    # kwargs for off_axis_projection
    oap_kwargs = {}
    oap_kwargs['weight'] = (None, 'cell_mass')
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
