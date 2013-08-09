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

use_tmpdir = False


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def setup_dir():
    if use_tmpdir:
        curdir = os.getcwd()
        # Perform I/O in safe place instead of yt main dir
        tmpdir = tempfile.mkdtemp()
        os.chdir(tmpdir)
    else:
        curdir, tmpdir = None, None
    return curdir, tmpdir


def teardown_dir(curdir, tmpdir):
    if use_tmpdir:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)


def setup_pf():
    # args for off_axis_projection
    test_pf = fake_random_pf(64)
    c = [0.5, 0.5, 0.5]
    norm = [0.5, 0.5, 0.5]
    W = test_pf.domain_width
    N = 64
    field = "Density"
    cam_args = [test_pf, c, norm, W, N, field]
    return cam_args


def setup_transfer_function(pf, camera_type):
    if camera_type in ['perspective', 'camera', 'stereopair', 'interactive']:
        mi, ma = pf.h.all_data().quantities['Extrema']('Density')[0]
        mi, ma = np.log10(mi), np.log10(ma)
        tf = ColorTransferFunction((mi-1., ma+1.), grey_opacity=True)
        Nsamples = 4
        tf.add_layers(Nsamples, w=0.02, col_bounds=(mi, ma),
                      alpha=np.logspace(1.0, 2.0, Nsamples), colormap='RdBu_r')
        #tf.map_to_colormap(mi,ma, scale=10., colormap='RdBu_r')
        return tf
    elif camera_type in ['healpix']:
        return ProjectionTransferFunction()
    else:
        pass


def test_camera():
    curdir, tmpdir = setup_dir()
    pf, c, L, W, N, field = setup_pf()
    tf = setup_transfer_function(pf, 'camera')

    cam = pf.h.camera(c, L, W, N, transfer_function=tf)
    cam.snapshot('camera.png')
    assert_fname('camera.png')
    teardown_dir(curdir, tmpdir)


def test_perspective_camera():
    curdir, tmpdir = setup_dir()
    pf, c, L, W, N, field = setup_pf()
    tf = setup_transfer_function(pf, 'camera')

    cam = PerspectiveCamera(c, L, W, N, pf=pf, transfer_function=tf)
    cam.snapshot('perspective.png')
    assert_fname('perspective.png')
    teardown_dir(curdir, tmpdir)


def test_interactive_camera():
    curdir, tmpdir = setup_dir()
    pf, c, L, W, N, field = setup_pf()
    tf = setup_transfer_function(pf, 'camera')

    cam = InteractiveCamera(c, L, W, N, pf=pf, transfer_function=tf)
    cam.snapshot('interactive.png')
    assert_fname('interactive.png')
    teardown_dir(curdir, tmpdir)


def test_projection_camera():
    curdir, tmpdir = setup_dir()
    pf, c, L, W, N, field = setup_pf()

    cam = ProjectionCamera(c, L, W, N, pf=pf, field='Density')
    cam.snapshot('projection.png')
    assert_fname('projection.png')
    teardown_dir(curdir, tmpdir)


def test_stereo_camera():
    curdir, tmpdir = setup_dir()
    pf, c, L, W, N, field = setup_pf()
    tf = setup_transfer_function(pf, 'camera')

    cam = pf.h.camera(c, L, W, N, transfer_function=tf)
    stereo_cam = StereoPairCamera(cam)
    # Take image
    cam1, cam2 = stereo_cam.split()
    cam1.snapshot(fn='stereo1.png')
    cam2.snapshot(fn='stereo2.png')
    assert_fname('stereo1.png')
    assert_fname('stereo2.png')
    teardown_dir(curdir, tmpdir)


def test_camera_movement():
    curdir, tmpdir = setup_dir()
    pf, c, L, W, N, field = setup_pf()
    tf = setup_transfer_function(pf, 'camera')

    cam = pf.h.camera(c, L, W, N, transfer_function=tf)
    cam.zoom(0.5)
    for snap in cam.zoomin(2.0, 3):
        snap
    for snap in cam.move_to(np.array(c) + 0.1, 3,
                            final_width=None, exponential=False):
        snap
    for snap in cam.move_to(np.array(c) - 0.1, 3,
                            final_width=2.0*W, exponential=False):
        snap
    for snap in cam.move_to(np.array(c), 3,
                            final_width=1.0*W, exponential=True):
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
    teardown_dir(curdir, tmpdir)
