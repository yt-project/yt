"""
Answer test to verify VR orientation and rotation is correct
"""

# -----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
import os
import tempfile

import numpy as np
import pytest

from yt.utilities.answer_testing.answer_tests import (
    VR_image_comparison,
    generic_image,
)
from yt.visualization.volume_rendering.api import off_axis_projection


@pytest.mark.answer_test
@pytest.mark.usefixtures("hashing")
class TestVROrientation:
    self.answer_file = None

    def test_vr_images(self, ds_vr, sc, lens_type):
        n_frames = 1
        theta = np.pi / n_frames
        cam = sc.add_camera(ds_vr, lens_type=lens_type)
        cam.resolution = (1000, 1000)
        cam.position = ds_vr.arr(np.array([-4.0, 0.0, 0.0]), "code_length")
        cam.switch_orientation(
            normal_vector=[1.0, 0.0, 0.0], north_vector=[0.0, 0.0, 1.0]
        )
        cam.set_width(ds_vr.domain_width * 2.0)
        test1 = VR_image_comparison(sc)
        self.hashes.update({"test1": test1})
        for i in range(n_frames):
            center = ds_vr.arr([0, 0, 0], "code_length")
            cam.yaw(theta, rot_center=center)
            test2 = VR_image_comparison(sc)
            # Updating nested dictionaries doesn't add the new key, it
            # overwrites the old one (so d.update({'key1' : {'subkey1' : 1}})
            # is d = {'key1' : {'subkey1' : 1}}. Then if you do
            # d.update({'key1' : {'subkey2' : 2}}), d = {'key1' : 'subkey2':2}},
            # so to add subkey2 to key1's subdictionary, you need to do
            # d['key1'].update({'subkey2' : 2}))
            if "test2" not in self.hashes:
                self.hashes.update({"test2": {str(i): test2}})
            else:
                self.hashes["test2"].update({str(i): test2})
        for i in range(n_frames):
            theta = np.pi / n_frames
            center = ds_vr.arr([0, 0, 0], "code_length")
            cam.pitch(theta, rot_center=center)
            test3 = VR_image_comparison(sc)
            if "test3" not in self.hashes:
                self.hashes.update({"test3": {str(i): test3}})
            else:
                self.hashes["test3"].update({str(i): test3})
        for i in range(n_frames):
            theta = np.pi / n_frames
            center = ds_vr.arr([0, 0, 0], "code_length")
            cam.roll(theta, rot_center=center)
            test4 = VR_image_comparison(sc)
            if "test4" not in self.hashes:
                self.hashes.update({"test4": {str(i): test4}})
            else:
                self.hashes["test4"].update({str(i): test4})

    def test_orientations(self, ds_vr, orientation):
        center = [0.5, 0.5, 0.5]
        width = [1.0, 1.0, 1.0]
        image = off_axis_projection(
            ds_vr, center, orientation, width, 512, "density", no_ghost=False
        )

        def offaxis_image_func():
            tmpfd, tmpfname = tempfile.mkstemp(suffix=".png")
            os.close(tmpfd)
            image.write_image(tmpfname)
            return tmpfname

        img_fname = offaxis_image_func()
        gi = generic_image(img_fname)
        self.hashes.update({"generic_image": gi})
