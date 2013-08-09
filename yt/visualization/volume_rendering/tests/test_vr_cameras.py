"""
Test for Volume Rendering Cameras, and their movement. 

Author: Samuel Skillman <samskillman@gmail.com>
Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Samuel Skillman.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import os
import os.path
import tempfile
import shutil
from yt.testing import \
    fake_random_pf
import numpy as np
from yt.mods import ColorTransferFunction, ProjectionTransferFunction
from yt.visualization.volume_rendering.api import \
    PerspectiveCamera, StereoPairCamera, InteractiveCamera, ProjectionCamera
from yt.visualization.tests.test_plotwindow import assert_fname
from unittest import TestCase

# This toggles using a temporary directory. Turn off to examine images.
use_tmpdir = True 


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


class CameraTest(TestCase):
    def setUp(self):
        if use_tmpdir:
            self.curdir = os.getcwd()
            # Perform I/O in safe place instead of yt main dir
            self.tmpdir = tempfile.mkdtemp()
            os.chdir(self.tmpdir)
        else:
            self.curdir, self.tmpdir = None, None

        self.pf, self.c, self.L, self.W, self.N, self.field = self.setup_pf()

    def tearDown(self):
        if use_tmpdir:
            os.chdir(self.curdir)
            shutil.rmtree(self.tmpdir)

    def setup_pf(self):
        # args for off_axis_projection
        test_pf = fake_random_pf(64)
        c = test_pf.domain_center
        norm = [0.5, 0.5, 0.5]
        W = 1.5*test_pf.domain_width
        N = 64
        field = "Density"
        cam_args = [test_pf, c, norm, W, N, field]
        return cam_args

    def setup_transfer_function(self, pf, camera_type):
        if camera_type in ['perspective', 'camera',
                           'stereopair', 'interactive']:
            mi, ma = pf.h.all_data().quantities['Extrema']('Density')[0]
            tf = ColorTransferFunction((mi-1., ma+1.), grey_opacity=True)
            tf.map_to_colormap(mi, ma, scale=10., colormap='RdBu_r')
            return tf
        elif camera_type in ['healpix']:
            return ProjectionTransferFunction()
        else:
            pass

    def test_camera(self):
        pf = self.pf
        tf = self.setup_transfer_function(pf, 'camera')
        cam = self.pf.h.camera(self.c, self.L, self.W, self.N,
                               transfer_function=tf)
        cam.snapshot('camera.png')
        assert_fname('camera.png')

    def test_source_camera(self):
        pf = self.pf
        tf = self.setup_transfer_function(pf, 'camera')
        source = pf.h.sphere(pf.domain_center, pf.domain_width[0]*0.5)

        cam = pf.h.camera(self.c, self.L, self.W, self.N,
                          transfer_function=tf, source=source)
        cam.snapshot('source_camera.png')
        assert_fname('source_camera.png')

    def test_perspective_camera(self):
        pf = self.pf
        tf = self.setup_transfer_function(pf, 'camera')

        cam = PerspectiveCamera(self.c, self.L, self.W, self.N, pf=pf,
                                transfer_function=tf)
        cam.snapshot('perspective.png')
        assert_fname('perspective.png')

    def test_interactive_camera(self):
        pf = self.pf
        tf = self.setup_transfer_function(pf, 'camera')

        cam = InteractiveCamera(self.c, self.L, self.W, self.N, pf=pf,
                                transfer_function=tf)
        # Can't take a snapshot here since IC uses pylab.'

    def test_projection_camera(self):
        pf = self.pf

        cam = ProjectionCamera(self.c, self.L, self.W, self.N, pf=pf,
                               field='Density')
        cam.snapshot('projection.png')
        assert_fname('projection.png')

    def test_stereo_camera(self):
        pf = self.pf
        tf = self.setup_transfer_function(pf, 'camera')

        cam = pf.h.camera(self.c, self.L, self.W, self.N, transfer_function=tf)
        stereo_cam = StereoPairCamera(cam)
        # Take image
        cam1, cam2 = stereo_cam.split()
        cam1.snapshot(fn='stereo1.png')
        cam2.snapshot(fn='stereo2.png')
        assert_fname('stereo1.png')
        assert_fname('stereo2.png')

    def test_camera_movement(self):
        pf = self.pf
        tf = self.setup_transfer_function(pf, 'camera')

        cam = pf.h.camera(self.c, self.L, self.W, self.N, transfer_function=tf)
        cam.zoom(0.5)
        for snap in cam.zoomin(2.0, 3):
            snap
        for snap in cam.move_to(np.array(self.c) + 0.1, 3,
                                final_width=None, exponential=False):
            snap
        for snap in cam.move_to(np.array(self.c) - 0.1, 3,
                                final_width=2.0*self.W, exponential=False):
            snap
        for snap in cam.move_to(np.array(self.c), 3,
                                final_width=1.0*self.W, exponential=True):
            snap
        cam.rotate(np.pi/10)
        cam.pitch(np.pi/10)
        cam.yaw(np.pi/10)
        cam.roll(np.pi/10)
        for snap in cam.rotation(np.pi, 3, rot_vector=None):
            snap
        for snap in cam.rotation(np.pi, 3, rot_vector=np.random.random(3)):
            snap
        cam.snapshot('final.png')
        assert_fname('final.png')
