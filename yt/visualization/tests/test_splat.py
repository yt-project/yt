"""
Test for write_bitmap and add_rgba_points



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
import numpy as np
import yt
from yt.testing import \
    assert_equal
from yt.utilities.lib.api import add_rgba_points_to_image


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_splat():
    """Tests functionality of off_axis_projection and write_projection."""
    # Perform I/O in safe place instead of yt main dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    N = 16 
    Np = int(1e2)
    image = np.zeros([N,N,4])
    xs = np.random.random(Np)
    ys = np.random.random(Np)

    cbx = yt.visualization.color_maps.mcm.RdBu
    cs = cbx(np.random.random(Np))
    add_rgba_points_to_image(image, xs, ys, cs)

    before_hash = image.copy()
    fn = 'tmp.png'
    yt.write_bitmap(image, fn)
    yield assert_equal, os.path.exists(fn), True
    os.remove(fn)
    yield assert_equal, before_hash, image

    os.chdir(curdir)
    # clean up
    shutil.rmtree(tmpdir)
