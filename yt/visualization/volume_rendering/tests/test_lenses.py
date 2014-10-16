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
from yt.mods import *
from yt.testing import \
    fake_random_ds
from yt.visualization.volume_rendering.scene import Scene
from yt.visualization.volume_rendering.camera import Camera
from yt.visualization.volume_rendering.render_source import VolumeSource
from time import time

#ds = fake_random_ds(8)
if False:
    ds = load('/home/skillman/kipac/data/IsolatedGalaxy/galaxy0030/galaxy0030')
    w = (ds.domain_width*30).in_units('code_length')
    w = ds.arr(w, 'code_length')
    sc = Scene()
    cam = Camera(ds, lens_type='perspective')
    print "WIDTH: ", w
    cam.set_width(w)
    vol = VolumeSource(ds, field=('gas', 'density'))
    sc.set_default_camera(cam)
    sc.add_source(vol)

    t = -time()
    sc.render('test_perspective.png', clip_ratio=None)
    t += time()
    print 'Total time: %e' % t

if True:
    ds = load('/home/skillman/kipac/data/IsolatedGalaxy/galaxy0030/galaxy0030')
    dd = ds.h.sphere(ds.domain_center, ds.domain_width[0] / 100)
    sc = Scene()
    cam = Camera(ds, lens_type='fisheye')
    cam.lens.fov = 360.0
    cam.resolution = (512,512)
    cam.set_width(ds.domain_width)
    v, c = ds.find_max('density')
    # c = ds.arr(c, 'code_length')
    p = ds.domain_center.copy()
    cam.set_position(c-0.0005*ds.domain_width)
    # ds.field_info[('gas','density')].take_log=False
    vol = VolumeSource(dd, field=('gas', 'density'))
    tf = vol.transfer_function
    tf.grey_opacity = True
    # tf.map_to_colormap(tf.x_bounds[0], tf.x_bounds[1], scale=3000.0, 
    #                    colormap='RdBu')
    sc.set_default_camera(cam)
    sc.add_source(vol)

    t = -time()
    sc.render('test_fisheye.png', clip_ratio=6.0)
    t += time()
    print 'Total time: %e' % t
