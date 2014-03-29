
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
pf = load('/home/skillman/kipac/data/IsolatedGalaxy/galaxy0030/galaxy0030')
ds = pf.h.sphere(pf.domain_center, pf.domain_width[0] / 50)


sc = Scene()
vol = VolumeSource(ds, field=('gas', 'density'))
cam = Camera(ds)
sc.set_default_camera(cam)
sc.add_source(vol)
sc.render('test.png')

im, sc2 = volume_render(ds, field=('gas', 'temperature'))
im.write_png('test2.png')
