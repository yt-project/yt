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
from collections import OrderedDict
import os
import tempfile

import numpy as np
import pytest

from yt.testing import fake_vr_orientation_test_ds
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils
from yt.visualization.volume_rendering.api import \
    Scene, \
    VolumeSource, \
    ColorTransferFunction, \
    off_axis_projection


# Answer file
answer_file = 'vr_orientation_answers.yaml'


@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestVROrientation(fw.AnswerTest):
    def test_orientation(self):
        ds = fake_vr_orientation_test_ds()
        sc = Scene()
        vol = VolumeSource(ds, field=('gas', 'density'))
        sc.add_source(vol)
        tf = vol.transfer_function
        tf = ColorTransferFunction((0.1, 1.0))
        tf.sample_colormap(1.0, 0.01, colormap="coolwarm")
        tf.sample_colormap(0.8, 0.01, colormap="coolwarm")
        tf.sample_colormap(0.6, 0.01, colormap="coolwarm")
        tf.sample_colormap(0.3, 0.01, colormap="coolwarm")
        n_frames = 1
        orientations = [[-0.3, -0.1, 0.8]]
        theta = np.pi / n_frames
        hd = OrderedDict()
        hd['vr_image_comparison'] = OrderedDict()
        hd['generic_image'] = OrderedDict()
        hd['vr_image_comparison']['test1'] = OrderedDict()
        hd['vr_image_comparison']['test2'] = OrderedDict()
        hd['vr_image_comparison']['test3'] = OrderedDict()
        hd['vr_image_comparison']['test4'] = OrderedDict()
        hd['vr_image_comparison']['test5'] = OrderedDict()
        for lens_type in ['plane-parallel', 'perspective']:
            frame = 0
            cam = sc.add_camera(ds, lens_type=lens_type)
            cam.resolution = (1000, 1000)
            cam.position = ds.arr(np.array([-4., 0., 0.]), 'code_length')
            cam.switch_orientation(normal_vector=[1., 0., 0.],
                                   north_vector=[0., 0., 1.])
            cam.set_width(ds.domain_width*2.)
            test1_hd = utils.generate_hash(self.VR_image_comparison_test(sc))
            hd['vr_image_comparison']['test1'][lens_type] = test1_hd
            hd['vr_image_comparison']['test2'][lens_type] = OrderedDict()
            hd['vr_image_comparison']['test3'][lens_type] = OrderedDict()
            hd['vr_image_comparison']['test4'][lens_type] = OrderedDict()
            for i in range(n_frames):
                frame += 1
                center = ds.arr([0, 0, 0], 'code_length')
                cam.yaw(theta, rot_center=center)
                test2_hd = utils.generate_hash(self.VR_image_comparison_test(sc))
                hd['vr_image_comparison']['test2'][lens_type][str(i)] = test2_hd
            for i in range(n_frames):
                frame += 1
                theta = np.pi / n_frames
                center = ds.arr([0, 0, 0], 'code_length')
                cam.pitch(theta, rot_center=center)
                test3_hd = utils.generate_hash(self.VR_image_comparison_test(sc))
                hd['vr_image_comparison']['test3'][lens_type][str(i)] = test3_hd
            for i in range(n_frames):
                frame += 1
                theta = np.pi / n_frames
                center = ds.arr([0, 0, 0], 'code_length')
                cam.roll(theta, rot_center=center)
                test4_hd = utils.generate_hash(self.VR_image_comparison_test(sc))
                hd['vr_image_comparison']['test4'][lens_type][str(i)] = test4_hd
        center = [0.5, 0.5, 0.5]
        width = [1.0, 1.0, 1.0]
        for i, orientation in enumerate(orientations):
            image = off_axis_projection(ds, center, orientation, width,
                                        512, "density", no_ghost=False)
            def offaxis_image_func():
                tmpfd, tmpfname = tempfile.mkstemp(suffix='.png')
                os.close(tmpfd)
                image.write_image(tmpfname)
                return tmpfname
            img_fname = offaxis_image_func()
            test5_hd = utils.generate_hash(self.generic_image_test(img_fname))
            hd['generic_image'][orientation] = test5_hd
        hd = {'orientation' : hd}
        utils.handle_hashes(self.save_dir, answer_file, hd, self.answer_store)
