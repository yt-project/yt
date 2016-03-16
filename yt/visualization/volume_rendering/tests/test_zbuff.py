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

import os
import tempfile
import shutil
from yt.testing import fake_random_ds
from yt.visualization.volume_rendering.api import \
    Scene, \
    ZBuffer, \
    VolumeSource, \
    OpaqueSource
from yt.testing import assert_almost_equal
import numpy as np
from unittest import TestCase


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


class ZBufferTest(TestCase):
    # This toggles using a temporary directory. Turn off to examine images.
    use_tmpdir = True

    def setUp(self):
        np.random.seed(0)
        if self.use_tmpdir:
            self.curdir = os.getcwd()
            # Perform I/O in safe place instead of yt main dir
            self.tmpdir = tempfile.mkdtemp()
            os.chdir(self.tmpdir)
        else:
            self.curdir, self.tmpdir = None, None

    def tearDown(self):
        if self.use_tmpdir:
            os.chdir(self.curdir)
            shutil.rmtree(self.tmpdir)


    def test_composite_vr(self):
        ds = fake_random_ds(64)
        dd = ds.sphere(ds.domain_center, 0.45*ds.domain_width[0])
        ds.field_info[ds.field_list[0]].take_log=False

        sc = Scene()
        cam = sc.add_camera(ds)
        cam.resolution = (512,512)
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

        im = sc.render()
        im.write_png("composite.png")
        return im

    def test_nonrectangular_add(self):
        rgba1 = np.ones((64, 1, 4))
        z1 = np.expand_dims(np.arange(64.), 1)

        rgba2 = np.zeros((64, 1, 4))
        z2 = np.expand_dims(np.arange(63., -1., -1.), 1)

        exact_rgba = np.concatenate((np.ones(32), np.zeros(32)))
        exact_rgba = np.expand_dims(exact_rgba, 1)
        exact_rgba = np.dstack((exact_rgba, exact_rgba, exact_rgba, exact_rgba))
        
        exact_z = np.concatenate((np.arange(32.), np.arange(31.,-1.,-1.)))
        exact_z = np.expand_dims(exact_z, 1)
        
        buff1 = ZBuffer(rgba1, z1)
        buff2 = ZBuffer(rgba2, z2)
        
        buff = buff1 + buff2
        
        assert_almost_equal(buff.rgba, exact_rgba)
        assert_almost_equal(buff.z, exact_z)

    def test_rectangular_add(self):
        rgba1 = np.ones((8, 8, 4))
        z1 = np.arange(64.)
        z1 = z1.reshape((8, 8))
        buff1 = ZBuffer(rgba1, z1)

        rgba2 = np.zeros((8, 8, 4))
        z2 = np.arange(63., -1., -1.)
        z2 = z2.reshape((8, 8))
        buff2 = ZBuffer(rgba2, z2)

        buff = buff1 + buff2

        exact_rgba = np.empty((8, 8, 4), dtype=np.float64)
        exact_rgba[0:4,0:8,:] = 1.0
        exact_rgba[4:8,0:8,:] = 0.0
        
        exact_z = np.concatenate((np.arange(32.), np.arange(31., -1., -1.)))
        exact_z = np.expand_dims(exact_z, 1)
        exact_z = exact_z.reshape(8, 8)

        assert_almost_equal(buff.rgba, exact_rgba)
        assert_almost_equal(buff.z, exact_z)
