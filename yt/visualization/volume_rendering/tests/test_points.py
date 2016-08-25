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
    VolumeSource, \
    PointSource
import numpy as np
from unittest import TestCase


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


class PointsVRTest(TestCase):
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

    def test_points_vr(self):
        ds = fake_random_ds(64)
        dd = ds.sphere(ds.domain_center, 0.45*ds.domain_width[0])
        ds.field_info[ds.field_list[0]].take_log=False

        sc = Scene()
        cam = sc.add_camera(ds)
        cam.resolution = (512,512)
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

        points_source = PointSource(vertices, colors=colors)
        sc.add_source(points_source)
        im = sc.render()
        im.write_png("points.png")
        return im
