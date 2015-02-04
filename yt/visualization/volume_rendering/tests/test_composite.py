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
    VolumeSource, OpaqueSource, LineSource, BoxSource
from yt.utilities.lib.misc_utilities import lines
from yt.data_objects.api import ImageArray
import numpy as np
np.random.seed(0)

def test_composite_vr():
    ds = fake_random_ds(64)
    dd = ds.sphere(ds.domain_center, 0.45*ds.domain_width[0])
    ds.field_info[ds.field_list[0]].take_log=False

    sc = Scene()
    cam = Camera(ds)
    cam.resolution = (512,512)
    sc.set_default_camera(cam)
    vr = VolumeSource(dd, field=ds.field_list[0])
    vr.transfer_function.clear()
    vr.transfer_function.grey_opacity=False
    vr.transfer_function.map_to_colormap(0.0, 1.0, scale=0.1, colormap="Reds")
    sc.add_source(vr)
    
    cam.set_width( 1.8*ds.domain_width )
    cam.lens.setup_box_properties(cam)

    # DRAW SOME LINES
    npoints = 100
    vertices = np.random.random([npoints, 3])
    colors = np.random.random([npoints, 4])
    colors[:,3] = 0.01

    box_source = BoxSource(ds.domain_left_edge, ds.domain_right_edge, color=[1.,1.,1.,0.01])
    sc.add_source(box_source)
    
    box_source = BoxSource(ds.domain_left_edge + np.array([0.1,0.,0.3])*ds.domain_left_edge.uq,
            ds.domain_right_edge-np.array([0.1,0.2,0.3])*ds.domain_left_edge.uq, 
            color=np.array([0.0, 1.0, 0.0, 0.01]))
    sc.add_source(box_source)
    
    line_source = LineSource(vertices, colors)
    sc.add_source(line_source)

    im = sc.composite()
    im = ImageArray(im.d)
    im.write_png("composite.png")
    return im

if __name__ == "__main__":
    im = test_composite_vr()
