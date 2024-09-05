import os
import tempfile

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pytest

from yt.testing import fake_vr_orientation_test_ds
from yt.visualization.volume_rendering.api import (
    Scene,
    create_volume_source,
    off_axis_projection,
)


def scene_to_mpl_figure(scene):
    """helper function to convert a scene image rendering to matplotlib
    so we can rely on pytest-mpl to compare images
    """
    tmpfd, tmpname = tempfile.mkstemp(suffix=".png")
    os.close(tmpfd)
    scene.save(tmpname, sigma_clip=1.0)
    image = mpl.image.imread(tmpname)
    os.remove(tmpname)

    fig, ax = plt.subplots()
    ax.set(aspect="equal")
    ax.imshow(image)
    return fig


class TestOrientation:
    @classmethod
    def setup_class(cls):
        cls.ds = fake_vr_orientation_test_ds()

        cls.scene = Scene()

        vol = create_volume_source(cls.ds, field=("gas", "density"))
        cls.scene.add_source(vol)

        cls._last_lense_type = None

    @classmethod
    def set_camera(cls, lens_type):
        # this method isn't thread-safe
        # if lens_type == cls._last_lense_type:
        #    return

        cls._last_lense_type = lens_type

        cam = cls.scene.add_camera(cls.ds, lens_type=lens_type)
        cam.resolution = (1000, 1000)
        cam.position = cls.ds.arr(np.array([-4.0, 0.0, 0.0]), "code_length")
        cam.switch_orientation(
            normal_vector=[1.0, 0.0, 0.0], north_vector=[0.0, 0.0, 1.0]
        )
        cam.set_width(cls.ds.domain_width * 2.0)
        cls.camera = cam

    @pytest.mark.parametrize("lens_type", ["perspective", "plane-parallel"])
    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_vr_orientation_lense_type(self, lens_type):
        # note that a previous version of this test proved flaky
        # and required a much lower precision for plane-parallel
        # https://github.com/yt-project/yt/issue/3069
        # https://github.com/yt-project/yt/pull/3068
        # https://github.com/yt-project/yt/pull/3294
        self.set_camera(lens_type)
        return scene_to_mpl_figure(self.scene)

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_vr_orientation_yaw(self):
        self.set_camera("plane-parallel")
        center = self.ds.arr([0, 0, 0], "code_length")
        self.camera.yaw(np.pi, rot_center=center)
        return scene_to_mpl_figure(self.scene)

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_vr_orientation_pitch(self):
        self.set_camera("plane-parallel")
        center = self.ds.arr([0, 0, 0], "code_length")
        self.camera.pitch(np.pi, rot_center=center)
        return scene_to_mpl_figure(self.scene)

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_vr_orientation_roll(self):
        self.set_camera("plane-parallel")
        center = self.ds.arr([0, 0, 0], "code_length")
        self.camera.roll(np.pi, rot_center=center)
        return scene_to_mpl_figure(self.scene)

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_vr_orientation_off_axis_projection(self):
        image = off_axis_projection(
            self.ds,
            center=[0.5, 0.5, 0.5],
            normal_vector=[-0.3, -0.1, 0.8],
            width=[1.0, 1.0, 1.0],
            resolution=512,
            item=("gas", "density"),
            no_ghost=False,
        )

        fig, ax = plt.subplots()
        ax.imshow(image)
        return fig
