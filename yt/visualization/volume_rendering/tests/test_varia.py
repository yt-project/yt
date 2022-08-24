import os
import shutil
import tempfile
from unittest import TestCase

import numpy as np

import yt
from yt.testing import fake_random_ds
from yt.visualization.volume_rendering.api import Scene, create_volume_source


def setup():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


class VariousVRTests(TestCase):
    # This toggles using a temporary directory. Turn off to examine images.
    use_tmpdir = True

    def setUp(self):
        if self.use_tmpdir:
            self.curdir = os.getcwd()
            # Perform I/O in safe place instead of yt main dir
            self.tmpdir = tempfile.mkdtemp()
            os.chdir(self.tmpdir)
        else:
            self.curdir, self.tmpdir = None, None

        self.ds = fake_random_ds(32)

    def tearDown(self):
        if self.use_tmpdir:
            os.chdir(self.curdir)
            shutil.rmtree(self.tmpdir)
        del self.ds

    def test_simple_scene_creation(self):
        yt.create_scene(self.ds)

    def test_modify_transfer_function(self):
        im, sc = yt.volume_render(self.ds)

        volume_source = sc.get_source(0)
        tf = volume_source.transfer_function
        tf.clear()
        tf.grey_opacity = True
        tf.add_layers(3, colormap="RdBu")
        sc.render()

    def test_multiple_fields(self):
        im, sc = yt.volume_render(self.ds)

        volume_source = sc.get_source(0)
        volume_source.set_field(("gas", "velocity_x"))
        volume_source.set_weight_field(("gas", "density"))
        sc.render()

    def test_rotation_volume_rendering(self):
        im, sc = yt.volume_render(self.ds)

        sc.camera.yaw(np.pi)
        sc.render()

    def test_simple_volume_rendering(self):
        im, sc = yt.volume_render(self.ds, sigma_clip=4.0)

    def test_lazy_volume_source_construction(self):
        sc = Scene()
        source = create_volume_source(self.ds.all_data(), ("gas", "density"))

        assert source._volume is None
        assert source._transfer_function is None

        source.tfh.bounds = (0.1, 1)

        source.set_log(False)

        assert not source.log_field
        assert source.transfer_function.x_bounds == [0.1, 1]
        assert source._volume is None

        source.set_log(True)

        assert source.log_field
        assert source.transfer_function.x_bounds == [-1, 0]
        assert source._volume is None

        source.transfer_function = None
        source.tfh.bounds = None

        ad = self.ds.all_data()

        np.testing.assert_allclose(
            source.transfer_function.x_bounds,
            np.log10(ad.quantities.extrema(("gas", "density"))),
        )
        assert source.tfh.log == source.log_field

        source.set_field(("gas", "velocity_x"))
        source.set_log(False)

        assert source.transfer_function.x_bounds == list(
            ad.quantities.extrema(("gas", "velocity_x"))
        )
        assert source._volume is None

        source.set_field(("gas", "density"))

        assert source.volume is not None
        assert not source.volume._initialized
        assert source.volume.fields is None

        del source.volume
        assert source._volume is None

        sc.add_source(source)

        sc.add_camera()

        sc.render()

        assert source.volume is not None
        assert source.volume._initialized
        assert source.volume.fields == [("gas", "density")]
        assert source.volume.log_fields == [True]

        source.set_field(("gas", "velocity_x"))
        source.set_log(False)

        sc.render()

        assert source.volume is not None
        assert source.volume._initialized
        assert source.volume.fields == [("gas", "velocity_x")]
        assert source.volume.log_fields == [False]
