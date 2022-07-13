import os
import shutil
import tempfile
import unittest
import warnings

import numpy as np

from yt.data_objects.image_array import ImageArray
from yt.testing import assert_equal, requires_module


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True
    np.seterr(all="ignore")


def dummy_image(kstep, nlayers):
    im = np.zeros([64, 128, nlayers])
    for i in range(im.shape[0]):
        for k in range(im.shape[2]):
            im[i, :, k] = np.linspace(0.0, kstep * k, im.shape[1])
    return im


def test_rgba_rescale():
    im_arr = ImageArray(dummy_image(10.0, 4))

    new_im = im_arr.rescale(inline=False)
    assert_equal(im_arr[:, :, :3].max(), 2 * 10.0)
    assert_equal(im_arr[:, :, 3].max(), 3 * 10.0)
    assert_equal(new_im[:, :, :3].sum(axis=2).max(), 1.0)
    assert_equal(new_im[:, :, 3].max(), 1.0)

    im_arr.rescale()
    assert_equal(im_arr[:, :, :3].sum(axis=2).max(), 1.0)
    assert_equal(im_arr[:, :, 3].max(), 1.0)

    im_arr.rescale(cmax=0.0, amax=0.0)
    assert_equal(im_arr[:, :, :3].sum(axis=2).max(), 1.0)
    assert_equal(im_arr[:, :, 3].max(), 1.0)


class TestImageArray(unittest.TestCase):

    tmpdir = None
    curdir = None

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def test_image_arry_units(self):
        im_arr = ImageArray(dummy_image(0.3, 3), units="cm")

        assert str(im_arr.units) == "cm"

        new_im = im_arr.in_units("km")

        assert str(new_im.units) == "km"

    @requires_module("h5py")
    def test_image_array_hdf5(self):
        myinfo = {
            "field": "dinosaurs",
            "east_vector": np.array([1.0, 0.0, 0.0]),
            "north_vector": np.array([0.0, 0.0, 1.0]),
            "normal_vector": np.array([0.0, 1.0, 0.0]),
            "width": 0.245,
            "type": "rendering",
        }

        im_arr = ImageArray(dummy_image(0.3, 3), units="cm", info=myinfo)
        im_arr.save("test_3d_ImageArray", png=False)

        im = np.zeros([64, 128])
        for i in range(im.shape[0]):
            im[i, :] = np.linspace(0.0, 0.3 * 2, im.shape[1])

        myinfo = {
            "field": "dinosaurs",
            "east_vector": np.array([1.0, 0.0, 0.0]),
            "north_vector": np.array([0.0, 0.0, 1.0]),
            "normal_vector": np.array([0.0, 1.0, 0.0]),
            "width": 0.245,
            "type": "rendering",
        }

        im_arr = ImageArray(im, info=myinfo, units="cm")
        im_arr.save("test_2d_ImageArray", png=False)

        im_arr.save("test_2d_ImageArray_ds", png=False, dataset_name="Random_DS")

    def test_image_array_rgb_png(self):
        im = np.zeros([64, 128])
        for i in range(im.shape[0]):
            im[i, :] = np.linspace(0.0, 0.3 * 2, im.shape[1])
        im_arr = ImageArray(im)
        im_arr.save("standard-image", hdf5=False)

        im_arr = ImageArray(dummy_image(10.0, 3))
        im_arr.save("standard-png", hdf5=False)

    def test_image_array_rgba_png(self):
        im_arr = ImageArray(dummy_image(10.0, 4))
        im_arr.write_png("standard")
        im_arr.write_png("non-scaled.png", rescale=False)
        im_arr.write_png("black_bg.png", background="black")
        im_arr.write_png("white_bg.png", background="white")
        im_arr.write_png("green_bg.png", background=[0.0, 1.0, 0.0, 1.0])
        im_arr.write_png("transparent_bg.png", background=None)
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            im_arr.write_png("clipped.png", clip_ratio=0.5)
            assert str(w[0].message).startswith(
                "The 'clip_ratio' keyword argument is a deprecated alias for 'sigma_clip'. "
                "Please use 'sigma_clip' directly."
            )

    def test_image_array_background(self):
        im_arr = ImageArray(dummy_image(10.0, 4))
        im_arr.rescale()
        new_im = im_arr.add_background_color([1.0, 0.0, 0.0, 1.0], inline=False)
        new_im.write_png("red_bg.png")
        im_arr.add_background_color("black")
        im_arr.write_png("black_bg2.png")

    def test_write_image(self):
        im_arr = ImageArray(dummy_image(10.0, 4))
        im_arr.write_image("with_cmap", cmap_name="hot")
        im_arr.write_image("channel_1.png", channel=1)

    def test_clipping_value(self):
        im_arr = ImageArray(dummy_image(10.0, 4))
        clip_val1 = im_arr._clipping_value(1)
        clip_val2 = im_arr._clipping_value(1, im=im_arr)
        assert clip_val2 == clip_val1

        clip_val3 = im_arr._clipping_value(6)
        assert clip_val3 > clip_val2

        im_arr[:] = 1.0  # std will be 0, mean will be 1, so clip value will be 1
        assert im_arr._clipping_value(1) == 1.0

    def tearDown(self):
        os.chdir(self.curdir)
        # clean up
        shutil.rmtree(self.tmpdir)
