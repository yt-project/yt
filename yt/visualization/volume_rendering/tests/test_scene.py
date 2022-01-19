import os
import shutil
import tempfile
from unittest import TestCase

import numpy as np

from yt.testing import assert_fname, fake_random_ds, fake_vr_orientation_test_ds
from yt.visualization.volume_rendering.api import (
    create_scene,
    create_volume_source,
    volume_render,
)


def setup():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


class RotationTest(TestCase):

    # This toggles using a temporary directory. Turn off to examine images.
    use_tmpdir = True

    def setUp(self):
        if self.use_tmpdir:
            self.curdir = os.getcwd()
            # Perform I/O in safe place instead of yt main dir
            self.tmpdir = tempfile.mkdtemp()
            os.chdir(self.tmpdir)
        else:
            self.curdir, self.tmpdir = None, None

    def tearDown(self):
        if self.use_tmpdir:
            os.chdir(self.curdir)
            shutil.rmtree(self.tmpdir)

    def test_rotation(self):
        ds = fake_random_ds(32)
        ds2 = fake_random_ds(32)
        dd = ds.sphere(ds.domain_center, ds.domain_width[0] / 2)
        dd2 = ds2.sphere(ds2.domain_center, ds2.domain_width[0] / 2)

        im, sc = volume_render(dd, field=("gas", "density"))
        im.write_png("test.png")

        vol = sc.get_source(0)
        tf = vol.transfer_function
        tf.clear()
        mi, ma = dd.quantities.extrema(("gas", "density"))
        mi = np.log10(mi)
        ma = np.log10(ma)
        mi_bound = ((ma - mi) * (0.10)) + mi
        ma_bound = ((ma - mi) * (0.90)) + mi
        tf.map_to_colormap(mi_bound, ma_bound, scale=0.01, colormap="Blues_r")

        vol2 = create_volume_source(dd2, field=("gas", "density"))
        sc.add_source(vol2)

        tf = vol2.transfer_function
        tf.clear()
        mi, ma = dd2.quantities.extrema(("gas", "density"))
        mi = np.log10(mi)
        ma = np.log10(ma)
        mi_bound = ((ma - mi) * (0.10)) + mi
        ma_bound = ((ma - mi) * (0.90)) + mi
        tf.map_to_colormap(mi_bound, ma_bound, scale=0.01, colormap="Reds_r")
        fname = "test_scene.pdf"
        sc.save(fname, sigma_clip=6.0)
        assert_fname(fname)

        fname = "test_rot.png"
        sc.camera.pitch(np.pi)
        sc.render()
        sc.save(fname, sigma_clip=6.0, render=False)
        assert_fname(fname)


def test_annotations():
    from matplotlib.image import imread

    curdir = os.getcwd()
    tmpdir = tempfile.mkdtemp()
    os.chdir(tmpdir)
    ds = fake_vr_orientation_test_ds(N=16)
    sc = create_scene(ds)
    sc.annotate_axes()
    sc.annotate_domain(ds)
    sc.render()
    # ensure that there are actually red, green, blue, and white pixels
    # in the image. see Issue #1595
    im = sc._last_render
    for c in ([1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1], [1, 1, 1, 1]):
        assert np.where((im == c).all(axis=-1))[0].shape[0] > 0
    sc[0].tfh.tf.add_layers(10, colormap="cubehelix")
    sc.save_annotated(
        "test_scene_annotated.png",
        text_annotate=[[(0.1, 1.05), "test_string"]],
    )
    image = imread("test_scene_annotated.png")
    assert image.shape == sc.camera.resolution + (4,)
    os.chdir(curdir)
    shutil.rmtree(tmpdir)
