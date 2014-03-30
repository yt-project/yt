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
    fake_random_pf
from yt.visualization.volume_rendering.scene import Scene
from yt.visualization.volume_rendering.camera import Camera
from yt.visualization.volume_rendering.render_source import VolumeSource

pf = fake_random_pf(64)
sc = Scene()
cam = Camera(pf, lens_type='perspective')
w = (pf.domain_width[0]*1000).in_units('code_length')
print "WIDTH: ", w
cam.set_width(w)
vol = VolumeSource(pf, field=('gas', 'density'))
vol.transfer_function.clear()
vol.transfer_function.grey_opacity = True
vol.transfer_function.map_to_colormap(-1., 0., scale=100000.)
sc.set_default_camera(cam)
sc.add_source(vol)
sc.render('test_perspective.png')
