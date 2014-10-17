
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
    fake_random_ds
from yt.visualization.volume_rendering.scene import Scene
from yt.visualization.volume_rendering.api import volume_render
from yt.visualization.volume_rendering.camera import Camera
from yt.visualization.volume_rendering.render_source import VolumeSource

#ds = fake_random_ds(64)
#ds = load('/home/skillman/kipac/data/IsolatedGalaxy/galaxy0030/galaxy0030')
ds = load('/home/skillman/kipac/data/enzo_cosmology_plus/DD0040/DD0040')
dd = ds.h.sphere(ds.domain_center, ds.domain_width[0] / 2)

im, sc = volume_render(dd, field=('gas', 'density'))
im.write_png('test.png')

sc.default_camera.resolution = (512, 512)
vol = sc.get_source(0)
tf = vol.transfer_function
tf.clear()
tf.map_to_colormap(-28.5, -28, scale=10.001, colormap='Blues_r')

ds = load('/home/skillman/kipac/data/enzo_cosmology_plus/DD0046/DD0046')
dd = ds.h.sphere(ds.domain_center, ds.domain_width[0] / 2)
#dd = ds.box(ds.domain_left_edge, ds.domain_right_edge/2.)
vol2 = VolumeSource(dd, field=('gas', 'density'))
sc.add_source(vol2)

tf = vol2.transfer_function
tf.clear()
#tf.map_to_colormap(6.0, 7.0, scale=0.01, colormap='Reds_r')
tf.map_to_colormap(-28.5, -28, scale=10.1, colormap='Reds_r')

sc.render('test_composite.png', clip_ratio=6.0)

nrot = 10
for i in range(10):
    sc.default_camera.pitch(2*np.pi/4/nrot)
    sc.render('test_rot_%04i.png' % i, clip_ratio=6.0)


