import os
import os.path
import shutil
import tempfile

import numpy as np

import yt
from yt.testing import assert_equal
from yt.utilities.lib.api import add_rgba_points_to_image


def setup():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def test_splat():
    # Perform I/O in safe place instead of yt main dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng = np.random.RandomState(0x4D3D3D3)
    N = 16
    Np = int(1e2)
    image = np.zeros([N, N, 4])
    xs = prng.random_sample(Np)
    ys = prng.random_sample(Np)

    cbx = yt.visualization.color_maps.mcm.RdBu
    cs = cbx(prng.random_sample(Np))
    add_rgba_points_to_image(image, xs, ys, cs)

    before_hash = image.copy()
    fn = "tmp.png"
    yt.write_bitmap(image, fn)
    assert_equal(os.path.exists(fn), True)
    os.remove(fn)
    assert_equal(before_hash, image)

    os.chdir(curdir)
    # clean up
    shutil.rmtree(tmpdir)
