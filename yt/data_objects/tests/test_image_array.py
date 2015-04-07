import numpy as np
import os
import tempfile
import shutil
import unittest
from yt.data_objects.image_array import ImageArray
from yt.testing import \
    assert_equal


def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"
    np.seterr(all='ignore')


def dummy_image(kstep, nlayers):
    im = np.zeros([64, 128, nlayers])
    for i in range(im.shape[0]):
        for k in range(im.shape[2]):
            im[i, :, k] = np.linspace(0.0, kstep * k, im.shape[1])
    return im


def test_rgba_rescale():
    im_arr = ImageArray(dummy_image(10.0, 4))

    new_im = im_arr.rescale(inline=False)
    yield assert_equal, im_arr[:, :, :3].max(), 2 * 10.
    yield assert_equal, im_arr[:, :, 3].max(), 3 * 10.
    yield assert_equal, new_im[:, :, :3].sum(axis=2).max(), 1.0
    yield assert_equal, new_im[:, :, 3].max(), 1.0

    im_arr.rescale()
    yield assert_equal, im_arr[:, :, :3].sum(axis=2).max(), 1.0
    yield assert_equal, im_arr[:, :, 3].max(), 1.0


class TestImageArray(unittest.TestCase):

    tmpdir = None
    curdir = None

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def test_image_arry_units(self):
        im_arr = ImageArray(dummy_image(0.3, 3), input_units='cm')

        assert str(im_arr.units) == 'cm'

        new_im = im_arr.in_units('km')

        assert str(new_im.units) == 'km'

    def test_image_array_hdf5(self):
        myinfo = {'field': 'dinosaurs', 'east_vector': np.array([1., 0., 0.]),
                  'north_vector': np.array([0., 0., 1.]),
                  'normal_vector': np.array([0., 1., 0.]),
                  'width': 0.245, 'type': 'rendering'}

        im_arr = ImageArray(dummy_image(0.3, 3), input_units='cm', info=myinfo)
        im_arr.save('test_3d_ImageArray')

        im = np.zeros([64, 128])
        for i in range(im.shape[0]):
            im[i, :] = np.linspace(0., 0.3 * 2, im.shape[1])

        myinfo = {'field': 'dinosaurs', 'east_vector': np.array([1., 0., 0.]),
                  'north_vector': np.array([0., 0., 1.]),
                  'normal_vector': np.array([0., 1., 0.]),
                  'width': 0.245, 'type': 'rendering'}

        im_arr = ImageArray(im, info=myinfo, input_units='cm')
        im_arr.save('test_2d_ImageArray')


    def test_image_array_rgb_png(self):
        im_arr = ImageArray(dummy_image(10.0, 3))
        im_arr.write_png('standard.png')

    def test_image_array_rgba_png(self):
        im_arr = ImageArray(dummy_image(10.0, 4))
        im_arr.write_png('standard.png')
        im_arr.write_png('non-scaled.png', rescale=False)
        im_arr.write_png('black_bg.png', background='black')
        im_arr.write_png('white_bg.png', background='white')
        im_arr.write_png('green_bg.png', background=[0., 1., 0., 1.])
        im_arr.write_png('transparent_bg.png', background=None)

    def test_image_array_background(self):
        im_arr = ImageArray(dummy_image(10.0, 4))
        im_arr.rescale()
        new_im = im_arr.add_background_color([1., 0., 0., 1.], inline=False)
        new_im.write_png('red_bg.png')
        im_arr.add_background_color('black')
        im_arr.write_png('black_bg2.png')

    def tearDown(self):
        os.chdir(self.curdir)
        # clean up
        shutil.rmtree(self.tmpdir)
