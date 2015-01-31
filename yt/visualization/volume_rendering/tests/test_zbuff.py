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
    vr.transfer_function.grey_opacity=True
    vr.transfer_function.map_to_colormap(0.0, 1.0, scale=10.0, colormap="Reds")
    sc.add_source(vr)
    
    cam.set_width( 1.8*ds.domain_width )
    cam.lens.setup_box_properties(cam)

    # Create Arbitrary Z-buffer
    empty = cam.lens.new_image(cam)
    z = np.empty(empty.shape[:2], dtype='float64')
    # Let's put a blue plane right through the center
    z[:] = cam.width[2] / 2.
    empty[:,:,2] = 1.0 # Set blue to 1's
    empty[:,:,3] = 1.0 # Set alpha to 1's
    zbuffer = ZBuffer(empty, z)
    zsource = OpaqueSource()
    zsource.set_zbuffer(zbuffer)
    sc.add_source(zsource)

    # Now we call sc.composite() instead of render.
    im = sc.composite()
    im = ImageArray(im.d)
    im.write_png("composite.png")
    return im

if __name__ == "__main__":
    im = test_composite_vr()
