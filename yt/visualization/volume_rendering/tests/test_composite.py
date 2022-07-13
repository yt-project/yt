import os
import shutil
import tempfile
from unittest import TestCase

import numpy as np

from yt.data_objects.api import ImageArray
from yt.testing import fake_random_ds
from yt.visualization.volume_rendering.api import (
    BoxSource,
    LineSource,
    Scene,
    create_volume_source,
)


def setup():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


class CompositeVRTest(TestCase):
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
        dd = ds.sphere(ds.domain_center, 0.45 * ds.domain_width[0])
        ds.field_info[ds.field_list[0]].take_log = False

        sc = Scene()
        cam = sc.add_camera(ds)
        cam.resolution = (512, 512)
        vr = create_volume_source(dd, field=ds.field_list[0])
        vr.transfer_function.clear()
        vr.transfer_function.grey_opacity = True
        vr.transfer_function.map_to_colormap(0.0, 1.0, scale=3.0, colormap="Reds")
        sc.add_source(vr)

        cam.set_width(1.8 * ds.domain_width)
        cam.lens.setup_box_properties(cam)

        # DRAW SOME LINES
        npoints = 100
        vertices = np.random.random([npoints, 2, 3])
        colors = np.random.random([npoints, 4])
        colors[:, 3] = 0.10

        box_source = BoxSource(
            ds.domain_left_edge, ds.domain_right_edge, color=[1.0, 1.0, 1.0, 1.0]
        )
        sc.add_source(box_source)

        LE = ds.domain_left_edge + np.array([0.1, 0.0, 0.3]) * ds.domain_left_edge.uq
        RE = ds.domain_right_edge - np.array([0.1, 0.2, 0.3]) * ds.domain_left_edge.uq
        color = np.array([0.0, 1.0, 0.0, 0.10])
        box_source = BoxSource(LE, RE, color=color)
        sc.add_source(box_source)

        line_source = LineSource(vertices, colors)
        sc.add_source(line_source)

        im = sc.render()
        im = ImageArray(im.d)
        im.write_png("composite.png")
