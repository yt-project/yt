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

import numpy as np

import yt
from yt.testing import \
    fake_random_ds

def test_simple_scene_creation():
    ds = fake_random_ds(32)
    yt.create_scene(ds)

def test_modify_transfer_function():
    ds = fake_random_ds(32)
    im, sc = yt.volume_render(ds)

    volume_source = sc.get_source(0)
    tf = volume_source.transfer_function
    tf.clear()
    tf.grey_opacity = True
    tf.add_layers(3, colormap='RdBu')
    sc.render()

def test_multiple_fields():
    ds = fake_random_ds(32)
    im, sc = yt.volume_render(ds)

    volume_source = sc.get_source(0)
    volume_source.set_field(('gas', 'velocity_x'))
    volume_source.build_default_transfer_function()
    sc.render()

def test_rotation_volume_rendering():
    ds = fake_random_ds(32)
    im, sc = yt.volume_render(ds)

    angle = 2*np.pi
    frames = 10
    for i in range(frames):
        sc.camera.yaw(angle/frames)
        sc.render()

def test_simple_volume_rendering():
    ds = fake_random_ds(32)
    im, sc = yt.volume_render(ds, sigma_clip=4.0)
