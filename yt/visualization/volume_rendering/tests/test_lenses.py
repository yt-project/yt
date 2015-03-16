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
import yt
from yt.testing import fake_random_ds
from yt.visualization.volume_rendering.api import Scene, Camera, VolumeSource
from time import time
field = ("gas", "density")

def test_prospective_lens():
    ds = fake_random_ds(32, fields = field)
    #ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    w = (ds.domain_width*30).in_units('code_length')
    w = ds.arr(w, 'code_length')
    sc = Scene()
    cam = Camera(ds, lens_type='perspective')
    cam.set_width(w)
    vol = VolumeSource(ds, field=field)
    tf = vol.transfer_function
    tf.grey_opacity = True
    sc.set_camera(cam)
    sc.add_source(vol)
    sc.render('test_perspective_%s.png' % field[1], clip_ratio=6.0)

def test_fisheye_lens():
    ds = fake_random_ds(32, fields = field)
    #ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    dd = ds.sphere(ds.domain_center, ds.domain_width[0] / 10)
    sc = Scene()
    cam = Camera(dd, lens_type='fisheye')
    cam.lens.fov = 360.0
    cam.set_width(ds.domain_width)
    v, c = ds.find_max('density')
    p = ds.domain_center.copy()
    cam.set_position(c-0.0005*ds.domain_width)
    vol = VolumeSource(dd, field=field)
    tf = vol.transfer_function
    tf.grey_opacity = True
    sc.set_camera(cam)
    sc.add_source(vol)
    sc.render('test_fisheye_%s.png' % field[1], clip_ratio=6.0)

def test_plane_lens():
    ds = fake_random_ds(32, fields = field)
    #ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    dd = ds.sphere(ds.domain_center, ds.domain_width[0] / 10)
    sc = Scene()
    cam = Camera(dd, lens_type='plane-parallel')
    cam.set_width(ds.domain_width*1e-2)
    v, c = ds.find_max('density')
    p = ds.domain_center.copy()
    vol = VolumeSource(dd, field=field)
    tf = vol.transfer_function
    tf.grey_opacity = True
    sc.set_camera(cam)
    sc.add_source(vol)
    sc.render('test_plane_%s.png' % field[1], clip_ratio=6.0)
