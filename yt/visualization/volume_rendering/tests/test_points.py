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
    VolumeSource, OpaqueSource, LineSource, BoxSource, PointsSource
from yt.utilities.lib.misc_utilities import lines
from yt.data_objects.api import ImageArray
import numpy as np
np.random.seed(0)

def test_points_vr():
    ds = fake_random_ds(64)
    dd = ds.sphere(ds.domain_center, 0.45*ds.domain_width[0])
    ds.field_info[ds.field_list[0]].take_log=False

    sc = Scene()
    cam = Camera(ds)
    cam.resolution = (512,512)
    sc.set_camera(cam)
    vr = VolumeSource(dd, field=ds.field_list[0])
    vr.transfer_function.clear()
    vr.transfer_function.grey_opacity=False
    vr.transfer_function.map_to_colormap(0.0, 1.0, scale=10., colormap="Reds")
    sc.add_source(vr)
    
    cam.set_width( 1.8*ds.domain_width )
    cam.lens.setup_box_properties(cam)

    # DRAW SOME POINTS
    npoints = 1000
    vertices = np.random.random([npoints, 3])
    colors = np.random.random([npoints, 4])
    colors[:,3] = 0.10

    points_source = PointsSource(vertices, colors=colors)
    sc.add_source(points_source)
    im = sc.render()
    im.write_png("points.png")
    return im

if __name__ == "__main__":
    im = test_points_vr()
