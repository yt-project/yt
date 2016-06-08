"""
Test for Volume Rendering Lenses.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import tempfile
import shutil
from yt.testing import fake_random_ds
from yt.visualization.volume_rendering.api import Scene, VolumeSource
import numpy as np
from unittest import TestCase

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


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
        self.ds = fake_random_ds(32, fields=self.field)
        self.ds.index

    def tearDown(self):
        if self.use_tmpdir:
            os.chdir(self.curdir)
            shutil.rmtree(self.tmpdir)

    def test_perspective_lens(self):
        sc = Scene()
        cam = sc.add_camera(self.ds, lens_type='perspective')
        cam.position = self.ds.arr(np.array([1.0, 1.0, 1.0]), 'code_length')
        vol = VolumeSource(self.ds, field=self.field)
        tf = vol.transfer_function
        tf.grey_opacity = True
        sc.add_source(vol)
        sc.render()
        sc.save('test_perspective_%s.png' % self.field[1], sigma_clip=6.0)

    def test_stereoperspective_lens(self):
        sc = Scene()
        cam = sc.add_camera(self.ds, lens_type='stereo-perspective')
        cam.resolution = [1024, 512]
        cam.position = self.ds.arr(np.array([0.7, 0.7, 0.7]), 'code_length')
        vol = VolumeSource(self.ds, field=self.field)
        tf = vol.transfer_function
        tf.grey_opacity = True
        sc.add_source(vol)
        sc.render()
        sc.save('test_stereoperspective_%s.png' % self.field[1], sigma_clip=6.0)

    def test_fisheye_lens(self):
        dd = self.ds.sphere(self.ds.domain_center,
                            self.ds.domain_width[0] / 10)
        sc = Scene()
        cam = sc.add_camera(dd, lens_type='fisheye')
        cam.lens.fov = 360.0
        cam.set_width(self.ds.domain_width)
        v, c = self.ds.find_max('density')
        cam.set_position(c-0.0005*self.ds.domain_width)
        vol = VolumeSource(dd, field=self.field)
        tf = vol.transfer_function
        tf.grey_opacity = True
        sc.add_source(vol)
        sc.render()
        sc.save('test_fisheye_%s.png' % self.field[1], sigma_clip=6.0)

    def test_plane_lens(self):
        dd = self.ds.sphere(self.ds.domain_center,
                            self.ds.domain_width[0] / 10)
        sc = Scene()
        cam = sc.add_camera(dd, lens_type='plane-parallel')
        cam.set_width(self.ds.domain_width*1e-2)
        v, c = self.ds.find_max('density')
        vol = VolumeSource(dd, field=self.field)
        tf = vol.transfer_function
        tf.grey_opacity = True
        sc.add_source(vol)
        sc.render()
        sc.save('test_plane_%s.png' % self.field[1], sigma_clip=6.0)

    def test_spherical_lens(self):
        sc = Scene()
        cam = sc.add_camera(self.ds, lens_type='spherical')
        cam.resolution = [512, 256]
        cam.position = self.ds.arr(np.array([0.6, 0.5, 0.5]), 'code_length')
        vol = VolumeSource(self.ds, field=self.field)
        tf = vol.transfer_function
        tf.grey_opacity = True
        sc.add_source(vol)
        sc.render()
        sc.save('test_spherical_%s.png' % self.field[1], sigma_clip=6.0)

    def test_stereospherical_lens(self):
        w = (self.ds.domain_width).in_units('code_length')
        w = self.ds.arr(w, 'code_length')
        sc = Scene()
        cam = sc.add_camera(self.ds, lens_type='stereo-spherical')
        cam.resolution = [512, 512]
        cam.position = self.ds.arr(np.array([0.6, 0.5, 0.5]), 'code_length')
        vol = VolumeSource(self.ds, field=self.field)
        tf = vol.transfer_function
        tf.grey_opacity = True
        sc.add_source(vol)
        sc.render()
        sc.save('test_stereospherical_%s.png' % self.field[1], sigma_clip=6.0)
