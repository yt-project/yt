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
from yt import testing
from yt.utilities.answer_testing.framework import \
    requires_answer_testing, \
    VRImageComparisonTest, \
    GenericImageTest
from yt.visualization.volume_rendering.api import \
    Scene, \
    VolumeSource, \
    ColorTransferFunction, \
    off_axis_projection


@requires_answer_testing()
def test_orientation():
    ds = testing.fake_vr_orientation_test_ds()

    sc = Scene()

    vol = VolumeSource(ds, field=('gas', 'density'))

    tf = vol.transfer_function
    tf = ColorTransferFunction((0.1, 1.0))
    tf.sample_colormap(1.0, 0.01, colormap="coolwarm")
    tf.sample_colormap(0.8, 0.01, colormap="coolwarm")
    tf.sample_colormap(0.6, 0.01, colormap="coolwarm")
    tf.sample_colormap(0.3, 0.01, colormap="coolwarm")

    n_frames = 5
    theta = np.pi / n_frames
    decimals = 12

    for lens_type in ['plane-parallel', 'perspective']:
        frame = 0

        cam = sc.add_camera(ds, lens_type=lens_type)
        cam.resolution = (1000, 1000)
        cam.position = ds.arr(np.array([-4., 0., 0.]), 'code_length')
        cam.switch_orientation(normal_vector=[1., 0., 0.],
                               north_vector=[0., 0., 1.])
        cam.set_width(ds.domain_width*2.)

        sc.add_source(vol)
        yield VRImageComparisonTest(
            sc, ds, '%s_%04d' % (lens_type, frame), decimals)

        for i in range(n_frames):
            frame += 1
            center = ds.arr([0, 0, 0], 'code_length')
            cam.yaw(theta, rot_center=center)
            yield VRImageComparisonTest(
                sc, ds, 'yaw_%s_%04d' % (lens_type, frame), decimals)

        for i in range(n_frames):
            frame += 1
            theta = np.pi / n_frames
            center = ds.arr([0, 0, 0], 'code_length')
            cam.pitch(theta, rot_center=center)
            yield VRImageComparisonTest(
                sc, ds, 'pitch_%s_%04d' % (lens_type, frame), decimals)

        for i in range(n_frames):
            frame += 1
            theta = np.pi / n_frames
            center = ds.arr([0, 0, 0], 'code_length')
            cam.roll(theta, rot_center=center)
            yield VRImageComparisonTest(
                sc, ds, 'roll_%s_%04d' % (lens_type, frame), decimals)

    orientations = [ [1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0],
                     [0.0, 0.0, 1.0],
                     [0.5, 0.4, 0.7],
                     [-0.3, -0.1, 0.8] ]
    center = [0.5, 0.5, 0.5]
    width = [1.0, 1.0, 1.0]

    for i, orientation in enumerate(orientations):
        image = off_axis_projection(ds, center, orientation, width,
                                    512, "density", no_ghost=False)

        def offaxis_image_func(filename_prefix):
            return image.write_image(filename_prefix)

        test = GenericImageTest(ds, offaxis_image_func, decimals)
        test.prefix = "oap_orientation_{}".format(i)
        yield test
