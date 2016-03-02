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

from yt import load_uniform_grid
from yt.utilities.answer_testing.framework import \
    requires_answer_testing, \
    VRImageComparisonTest, \
    GenericImageTest
from yt.visualization.volume_rendering.api import \
    Scene, \
    Camera, \
    VolumeSource, \
    ColorTransferFunction, \
    off_axis_projection


def setup_ds():

    N = 96

    xmin = ymin = zmin = -1.0
    xmax = ymax = zmax = 1.0

    dcoord = (xmax - xmin)/N

    arr = np.zeros((N, N, N), dtype=np.float64)
    arr[:, :, :] = 1.e-4

    bbox = np.array([[xmin, xmax], [ymin, ymax], [zmin, zmax]])

    # coordinates -- in the notation data[i, j, k]
    x = (np.arange(N) + 0.5)*(xmax - xmin)/N + xmin
    y = (np.arange(N) + 0.5)*(ymax - ymin)/N + ymin
    z = (np.arange(N) + 0.5)*(zmax - zmin)/N + zmin

    x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")

    # sphere at the origin
    c = np.array([0.5*(xmin + xmax), 0.5*(ymin + ymax), 0.5*(zmin + zmax)])

    r = np.sqrt((x3d - c[0])**2 + (y3d - c[1])**2 + (z3d - c[2])**2)
    arr[r < 0.05] = 1.0

    arr[abs(x3d - xmin) < 2*dcoord] = 0.3
    arr[abs(y3d - ymin) < 2*dcoord] = 0.3
    arr[abs(z3d - zmin) < 2*dcoord] = 0.3

    # single cube on +x
    xc = 0.75
    dx = 0.05
    idx = np.logical_and(np.logical_and(x3d > xc-dx, x3d < xc+dx),
                         np.logical_and(np.logical_and(y3d > -dx, y3d < dx),
                                        np.logical_and(z3d > -dx, z3d < dx)))

    arr[idx] = 1.0

    # two cubes on +y
    dy = 0.05
    for yc in [0.65, 0.85]:

        idx = np.logical_and(np.logical_and(y3d > yc-dy, y3d < yc+dy),
                             np.logical_and(np.logical_and(x3d > -dy, x3d < dy),
                                            np.logical_and(z3d > -dy, z3d < dy)))

        arr[idx] = 0.8

    # three cubes on +z
    dz = 0.05
    for zc in [0.5, 0.7, 0.9]:

        idx = np.logical_and(np.logical_and(z3d > zc-dz, z3d < zc+dz),
                             np.logical_and(np.logical_and(x3d > -dz, x3d < dz),
                                            np.logical_and(y3d > -dz, y3d < dz)))

        arr[idx] = 0.6

    data = dict(density=(arr, "g/cm**3"))
    ds = load_uniform_grid(data, arr.shape, bbox=bbox)

    return ds


@requires_answer_testing()
def test_orientation():
    ds = setup_ds()

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

        cam = Camera(ds, lens_type='plane-parallel')
        cam.resolution = (1000, 1000)
        cam.position = ds.arr(np.array([-4., 0., 0.]), 'code_length')
        cam.switch_orientation(normal_vector=[1., 0., 0.],
                               north_vector=[0., 0., 1.])
        cam.set_width(ds.domain_width*2.)

        sc.camera = cam
        sc.add_source(vol)
        yield VRImageComparisonTest(
            sc, ds, '%s_%04d' % (lens_type, frame), decimals)

        for i in range(n_frames):
            frame += 1
            center = ds.arr([0, 0, 0], 'code_length')
            cam.yaw(theta, rot_center=center)
            sc.camera = cam
            yield VRImageComparisonTest(
                sc, ds, 'yaw_%s_%04d' % (lens_type, frame), decimals)

        for i in range(n_frames):
            frame += 1
            theta = np.pi / n_frames
            center = ds.arr([0, 0, 0], 'code_length')
            cam.pitch(theta, rot_center=center)
            sc.camera = cam
            yield VRImageComparisonTest(
                sc, ds, 'pitch_%s_%04d' % (lens_type, frame), decimals)

        for i in range(n_frames):
            frame += 1
            theta = np.pi / n_frames
            center = ds.arr([0, 0, 0], 'code_length')
            cam.roll(theta, rot_center=center)
            sc.camera = cam
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
