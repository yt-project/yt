import os
import shutil
import tempfile
from unittest import TestCase

import numpy as np

from yt.testing import fake_random_ds
from yt.visualization.volume_rendering.api import Scene, create_volume_source


def setup():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


class LensTest(TestCase):
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

        self.field = ("gas", "density")
        self.ds = fake_random_ds(32, fields=(self.field,), units=("g/cm**3",))
        self.ds.index

    def tearDown(self):
        if self.use_tmpdir:
            os.chdir(self.curdir)
            shutil.rmtree(self.tmpdir)

    def test_perspective_lens(self):
        sc = Scene()
        cam = sc.add_camera(self.ds, lens_type="perspective")
        cam.position = self.ds.arr(np.array([1.0, 1.0, 1.0]), "code_length")
        vol = create_volume_source(self.ds, field=self.field)
        tf = vol.transfer_function
        tf.grey_opacity = True
        sc.add_source(vol)
        sc.save(f"test_perspective_{self.field[1]}.png", sigma_clip=6.0)

    def test_stereoperspective_lens(self):
        sc = Scene()
        cam = sc.add_camera(self.ds, lens_type="stereo-perspective")
        cam.resolution = [256, 128]
        cam.position = self.ds.arr(np.array([0.7, 0.7, 0.7]), "code_length")
        vol = create_volume_source(self.ds, field=self.field)
        tf = vol.transfer_function
        tf.grey_opacity = True
        sc.add_source(vol)
        sc.save(f"test_stereoperspective_{self.field[1]}.png", sigma_clip=6.0)

    def test_fisheye_lens(self):
        dd = self.ds.sphere(self.ds.domain_center, self.ds.domain_width[0] / 10)
        sc = Scene()
        cam = sc.add_camera(dd, lens_type="fisheye")
        cam.lens.fov = 360.0
        cam.set_width(self.ds.domain_width)
        v, c = self.ds.find_max(("gas", "density"))
        cam.set_position(c - 0.0005 * self.ds.domain_width)
        vol = create_volume_source(dd, field=self.field)
        tf = vol.transfer_function
        tf.grey_opacity = True
        sc.add_source(vol)
        sc.save(f"test_fisheye_{self.field[1]}.png", sigma_clip=6.0)

    def test_plane_lens(self):
        dd = self.ds.sphere(self.ds.domain_center, self.ds.domain_width[0] / 10)
        sc = Scene()
        cam = sc.add_camera(dd, lens_type="plane-parallel")
        cam.set_width(self.ds.domain_width * 1e-2)
        v, c = self.ds.find_max(("gas", "density"))
        vol = create_volume_source(dd, field=self.field)
        tf = vol.transfer_function
        tf.grey_opacity = True
        sc.add_source(vol)
        sc.save(f"test_plane_{self.field[1]}.png", sigma_clip=6.0)

    def test_spherical_lens(self):
        sc = Scene()
        cam = sc.add_camera(self.ds, lens_type="spherical")
        cam.resolution = [256, 128]
        cam.position = self.ds.arr(np.array([0.6, 0.5, 0.5]), "code_length")
        vol = create_volume_source(self.ds, field=self.field)
        tf = vol.transfer_function
        tf.grey_opacity = True
        sc.add_source(vol)
        sc.save(f"test_spherical_{self.field[1]}.png", sigma_clip=6.0)

    def test_stereospherical_lens(self):
        w = (self.ds.domain_width).in_units("code_length")
        w = self.ds.arr(w, "code_length")
        sc = Scene()
        cam = sc.add_camera(self.ds, lens_type="stereo-spherical")
        cam.resolution = [256, 256]
        cam.position = self.ds.arr(np.array([0.6, 0.5, 0.5]), "code_length")
        vol = create_volume_source(self.ds, field=self.field)
        tf = vol.transfer_function
        tf.grey_opacity = True
        sc.add_source(vol)
        sc.save(f"test_stereospherical_{self.field[1]}.png", sigma_clip=6.0)
