
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
    create_volume_rendering
from yt.visualization.volume_rendering.lens import PlaneParallelLens
from yt.visualization.volume_rendering.camera import Camera
from yt.visualization.volume_rendering.render_source import VolumeSource

#pf = fake_random_pf(64)
pf = load('/home/skillman/kipac/data/IsolatedGalaxy/galaxy0030/galaxy0030')
ds = pf.h.sphere(pf.domain_center, pf.domain_width[0] / 50)


sc = Scene()
vol = VolumeSource(ds, field=('gas', 'density'))
cam = Camera()
cam.resolution = (512, 512)
lens = PlaneParallelLens()
sc.camera = cam

vol.build_defaults()
vol.transfer_function.grey_opacity=False
sc.add_source(vol)
cam.set_defaults_from_data_source(ds)
cam.set_lens(lens)

lens.expose(sc, cam, vol)
vol.current_image.write_png('test.png', clip_ratio=6.0)

vol.set_field(('io','Density'))
vol.build_defaults()
lens.expose(sc, cam, vol)
vol.current_image.write_png('test_op.png', clip_ratio=6.0)

#sc.render()


#sc = create_volume_rendering(ds, field=('gas', 'density'))
#sc.render('test.png')
#
#h = sc.get_handle()
#h.source.transfer_function.grey_opacity = True
#h.source.transfer_function.map_to_colormap(-2, 0.0, scale=50.0, colormap='RdBu_r')
#
#cam = h.camera
#for i in range(36):
#    cam.pitch(-2*np.pi / 36.)
#    sc.render('test_%04i.png' % i)

