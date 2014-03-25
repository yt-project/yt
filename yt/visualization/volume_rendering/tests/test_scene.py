
"""
Test for Volume Rendering Scene, and their movement.




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
    fake_random_pf
from yt.visualization.volume_rendering.scene import Scene, RenderScene, \
    create_volume_rendering
from yt.visualization.volume_rendering.camera import Camera
from yt.visualization.volume_rendering.render_source import VolumeSource


# def test_default(pf):
#     sc = RenderScene(pf)
#     for k, im in sc.render().iteritems():
#         write_bitmap(im, 'scene_%s.png' % k)
#     return sc
# 
# 
# def test_data_source(pf):
#     sc = Scene()
#     ds = pf.h.sphere(pf.domain_center, pf.domain_width[0] / 3)
#     vol = VolumeSource(ds, field=('gas', 'density'))
#     cam = Camera(ds)
#     sc.set_camera(cam)
#     sc.add_source(vol)
# 
#     vol.build_defaults()
#     for k, im in sc.render().iteritems():
#         write_bitmap(im, 'data_scene_%s.png' % k)
# 
# 
# def test_two_source(pf):
#     sc = Scene()
#     vol = VolumeSource(pf.h.sphere(pf.domain_center, pf.domain_width[0] / 3),
#                        field=('gas', 'density'))
#     sc.add_source(vol)
#     vol.build_defaults()
# 
#     vol = VolumeSource(pf.h.sphere(pf.domain_center / 3,
#                                    pf.domain_width[0] / 3),
#                        field=('gas', 'density'))
#     sc.add_source(vol)
#     for k, im in sc.render().iteritems():
#         write_bitmap(im, 'muliple_scene_%s.png' % k)
# 

#pf = load('/home/skillman/kipac/data/IsolatedGalaxy/galaxy0030/galaxy0030')
pf = fake_random_pf(64)
ds = pf.h.sphere(pf.domain_center, pf.domain_width[0] / 2)
sc = create_volume_rendering(ds, field=('gas', 'density'))
sc.render('test.png')
h = sc.get_handle()
h.source.transfer_function.grey_opacity = True
h.source.transfer_function.map_to_colormap(-0.5, -0.2, scale=50.0, colormap='RdBu_r')
cam = h.camera
for i in range(10):
    cam.yaw(np.pi / 10.)
    sc.render('test_%04i.png' % i)

