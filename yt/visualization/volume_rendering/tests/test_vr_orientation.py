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


import numpy as np
from nose.plugins.attrib import attr
from yt.testing import ANSWER_TEST_TAG, fake_vr_orientation_test_ds
from yt.utilities.answer_testing.framework import \
    VRImageComparisonTest, \
    GenericImageTest
from yt.visualization.volume_rendering.api import \
    Scene, \
    VolumeSource, \
    ColorTransferFunction, \
    off_axis_projection


@attr(ANSWER_TEST_TAG)
def test_orientation():
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
    decimals = 12
    test_name = "vr_orientation"

    for lens_type in ['plane-parallel', 'perspective']:
        frame = 0

        cam = sc.add_camera(ds, lens_type=lens_type)
        cam.resolution = (1000, 1000)
        cam.position = ds.arr(np.array([-4., 0., 0.]), 'code_length')
        cam.switch_orientation(normal_vector=[1., 0., 0.],
                               north_vector=[0., 0., 1.])
        cam.set_width(ds.domain_width*2.)
        desc = '%s_%04d' % (lens_type, frame)
        test1 = VRImageComparisonTest(sc, ds, desc, decimals)
        test1.answer_name = test_name
        yield test1

        for i in range(n_frames):
            frame += 1
            center = ds.arr([0, 0, 0], 'code_length')
            cam.yaw(theta, rot_center=center)
            desc = 'yaw_%s_%04d' % (lens_type, frame)
            test2 = VRImageComparisonTest(sc, ds, desc, decimals)
            test2.answer_name = test_name
            yield test2

        for i in range(n_frames):
            frame += 1
            theta = np.pi / n_frames
            center = ds.arr([0, 0, 0], 'code_length')
            cam.pitch(theta, rot_center=center)
            desc = 'pitch_%s_%04d' % (lens_type, frame)
            test3 = VRImageComparisonTest(sc, ds, desc, decimals)
            test3.answer_name = test_name
            yield test3

        for i in range(n_frames):
            frame += 1
            theta = np.pi / n_frames
            center = ds.arr([0, 0, 0], 'code_length')
            cam.roll(theta, rot_center=center)
            desc = 'roll_%s_%04d' % (lens_type, frame)
            test4 = VRImageComparisonTest(sc, ds, desc, decimals)
            test4.answer_name = test_name
            yield test4

    center = [0.5, 0.5, 0.5]
    width = [1.0, 1.0, 1.0]

    for i, orientation in enumerate(orientations):
        image = off_axis_projection(ds, center, orientation, width,
                                    512, "density", no_ghost=False)

        def offaxis_image_func(filename_prefix):
            return image.write_image(filename_prefix)

        test5 = GenericImageTest(ds, offaxis_image_func, decimals)
        test5.prefix = "oap_orientation_{}".format(i)
        test5.answer_name = test_name
        yield test5
