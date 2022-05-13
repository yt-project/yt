import os
import shutil
import tempfile
import unittest

import numpy as np
from nose.tools import assert_raises

from yt.testing import assert_equal, fake_random_ds
from yt.visualization.image_writer import (
    apply_colormap,
    multi_image_composite,
    splat_points,
    write_bitmap,
)


class TestImageWriter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.curdir = os.getcwd()
        os.chdir(cls.tmpdir)

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.curdir)
        shutil.rmtree(cls.tmpdir)

    def test_multi_image_composite(self):
        ds = fake_random_ds(64, nprocs=4, particles=16**3)
        center = [0.5, 0.5, 0.5]
        normal = [1, 1, 1]
        cut = ds.cutting(normal, center)
        frb = cut.to_frb((0.75, "unitary"), 64)
        multi_image_composite(
            "multi_channel1.png", frb[("index", "x")], frb[("index", "y")]
        )

        # Test multi_image_composite with user specified scaling values
        mi = ds.quan(0.1, "code_length")
        ma = ds.quan(0.9, "code_length")
        multi_image_composite(
            "multi_channel2.png",
            (frb[("index", "x")], mi, ma),
            [frb[("index", "y")], mi, None],
            green_channel=frb[("index", "z")],
            alpha_channel=frb[("gas", "density")],
        )

        # Test with numpy integer array
        x = np.array(np.random.randint(0, 256, size=(10, 10)), dtype="uint8")
        y = np.array(np.random.randint(0, 256, size=(10, 10)), dtype="uint8")
        multi_image_composite("multi_channel3.png", x, y)

    def test_write_bitmap(self):
        image = np.zeros([16, 16, 4], dtype="uint8")
        xs = np.random.rand(100)
        ys = np.random.rand(100)
        image = splat_points(image, xs, ys)
        png_str = write_bitmap(image.copy(), None)

        image_trans = image.swapaxes(0, 1).copy(order="C")
        png_str_trans = write_bitmap(image_trans, None, transpose=True)
        assert_equal(png_str, png_str_trans)

        with assert_raises(RuntimeError) as ex:
            write_bitmap(np.ones([16, 16]), None)
        desired = "Expecting image array of shape (N,M,3) or (N,M,4), received (16, 16)"
        assert_equal(str(ex.exception)[:50], desired[:50])

    def test_apply_colormap(self):
        x = np.array(np.random.randint(0, 256, size=(10, 10)), dtype="uint8")
        apply_colormap(x, color_bounds=None, cmap_name=None, func=lambda x: x**2)
