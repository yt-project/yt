
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
from yt.visualization.volume_rendering.scene import Scene, \
    volume_render
from yt.visualization.volume_rendering.camera import Camera
from yt.visualization.volume_rendering.render_source import VolumeSource

#pf = fake_random_pf(64)
#pf = load('/home/skillman/kipac/data/IsolatedGalaxy/galaxy0030/galaxy0030')
pf = load('/home/skillman/kipac/data/enzo_cosmology_plus/DD0046/DD0046')
ds = pf.h.sphere(pf.domain_center, pf.domain_width[0] / 2)

im, sc = volume_render(ds, field=('gas', 'density'))
im.write_png('test.png')

sc.default_camera.resolution = (512, 512)
vol = sc.get_source(0)
tf = vol.transfer_function
tf.clear()
tf.map_to_colormap(-30.5, -28, scale=0.001, colormap='Blues_r')

vol2 = VolumeSource(ds, field=('gas', 'temperature'))
sc.add_source(vol2)

tf = vol2.transfer_function
tf.clear()
tf.map_to_colormap(6.0, 7.0, scale=0.01, colormap='Reds_r')

sc.render('test_composite.png', clip_ratio=6.0)

nrot = 10
for i in range(10):
    sc.default_camera.pitch(2*np.pi/4/nrot)
    sc.render('test_rot_%04i.png' % i, clip_ratio=6.0)


