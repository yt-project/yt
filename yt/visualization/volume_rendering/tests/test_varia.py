"""
Miscellaneous tests for VR infrastructure
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import tempfile
import shutil
import numpy as np

import yt
from yt.testing import \
    fake_random_ds
from unittest import TestCase

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

class VariousVRTests(TestCase):
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

        self.ds = fake_random_ds(32)
    
    def tearDown(self):
        if self.use_tmpdir:
            os.chdir(self.curdir)
            shutil.rmtree(self.tmpdir)

    def test_simple_scene_creation(self):
        yt.create_scene(self.ds)

    def test_modify_transfer_function(self):
        im, sc = yt.volume_render(self.ds)

        volume_source = sc.get_source(0)
        tf = volume_source.transfer_function
        tf.clear()
        tf.grey_opacity = True
        tf.add_layers(3, colormap='RdBu')
        sc.render()

    def test_multiple_fields(self):
        im, sc = yt.volume_render(self.ds)

        volume_source = sc.get_source(0)
        volume_source.set_field(('gas', 'velocity_x'))
        volume_source.set_weight_field(('gas', 'density'))
        sc.render()

    def test_rotation_volume_rendering(self):
        im, sc = yt.volume_render(self.ds)

        angle = 2 * np.pi
        frames = 10
        for i in range(frames):
            sc.camera.yaw(angle / frames)
            sc.render()

    def test_simple_volume_rendering(self):
        im, sc = yt.volume_render(self.ds, sigma_clip=4.0)
