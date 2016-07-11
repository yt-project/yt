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
from yt.visualization.volume_rendering.render_source import VolumeSource
from yt.visualization.volume_rendering.scene import Scene
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

    def test_lazy_volume_source_construction(self):
        sc = Scene()
        source = VolumeSource(self.ds.all_data(), 'density')

        assert source._volume is None
        assert source._transfer_function is None

        source.tfh.bounds = (0.1, 1)

        source.set_log(False)

        assert source.log_field is False
        assert source.transfer_function.x_bounds == [0.1, 1]
        assert source._volume is None

        source.set_log(True)

        assert source.log_field is True
        assert source.transfer_function.x_bounds == [-1, 0]
        assert source._volume is None

        source.transfer_function = None
        source.tfh.bounds = None

        ad = self.ds.all_data()

        assert source.transfer_function.x_bounds == \
            list(np.log10(ad.quantities.extrema('density')))
        assert source.tfh.log == source.log_field

        source.set_field('velocity_x')
        source.set_log(False)

        assert source.transfer_function.x_bounds == \
            list(ad.quantities.extrema('velocity_x'))
        assert source._volume is None

        source.set_field('density')

        assert source.volume is not None
        assert source.volume._initialized is False
        assert source.volume.fields is None

        del source.volume
        assert source._volume is None

        sc.add_source(source)

        sc.add_camera()

        sc.render()

        assert source.volume is not None
        assert source.volume._initialized is True
        assert source.volume.fields == [('gas', 'density')]
        assert source.volume.log_fields == [True]

        source.set_field('velocity_x')
        source.set_log(False)

        sc.render()

        assert source.volume is not None
        assert source.volume._initialized is True
        assert source.volume.fields == [('gas', 'velocity_x')]
        assert source.volume.log_fields == [False]
