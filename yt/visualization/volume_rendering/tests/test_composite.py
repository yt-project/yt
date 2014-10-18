"""
Test for Composite VR.
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
from yt.visualization.volume_rendering.api import Scene, Camera, ZBuffer, \
    VolumeSource, OpaqueSource
from yt.utilities.lib.misc_utilities import lines
from yt.data_objects.api import ImageArray
import numpy as np
np.random.seed(0)

def test_composite_vr():
    ds = fake_random_ds(64)
    dd = ds.sphere(ds.domain_center, ds.domain_width[0] / 3.)
    ds.field_info[ds.field_list[0]].take_log=False

    sc = Scene()
    cam = Camera(ds)
    sc.set_default_camera(cam)
    vr = VolumeSource(dd, field=ds.field_list[0])
    vr.transfer_function.clear()
    vr.transfer_function.grey_opacity=True
    vr.transfer_function.map_to_colormap(0.0, 1.0, scale=3.0, colormap="Reds")
    sc.add_source(vr)

    op = OpaqueSource()
    empty = 0.0 * sc.default_camera.lens.new_image(cam)
    z = np.ones(empty.shape[:2]) * np.inf
    
    zbuff = ZBuffer(empty, z)
    op.set_zbuffer(zbuff)
    
    sc.add_source(op)
    cam.set_width( ds.domain_width )
    
    # DRAW SOME LINES
    npoints = 100
    vertices = 0.5 * np.random.random([npoints, 3])
    cam.lens.setup_box_properties(cam)
    px, py, dz = cam.lens.project_to_plane(cam, vertices)
    colors = np.random.random([npoints, 4])
    colors[:,3] = 1.0 
    lines(empty, px, py, colors, 24)
    dummy = -np.ones_like(empty)
    lines(dummy, px, py, np.vstack([dz, dz]), 24)
    z[:,:] = dummy[:,:,3]
    z[z==-1] = np.inf
    
    im = sc.composite()
    im = ImageArray(im.d)
    im.write_png("composite.png")
