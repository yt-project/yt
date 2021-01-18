import numpy as np

import yt.units as u
from yt.testing import assert_almost_equal, assert_equal, fake_random_ds
from yt.visualization.volume_rendering.api import Scene

valid_lens_types = [
    "plane-parallel",
    "perspective",
    "stereo-perspective",
    "fisheye",
    "spherical",
    "stereo-spherical",
]


def test_scene_and_camera_attributes():
    ds = fake_random_ds(64, length_unit=2, bbox=np.array([[-1, 1], [-1, 1], [-1, 1]]))
    sc = Scene()
    cam = sc.add_camera(ds)

    # test that initial values are correct in code units
    assert_equal(cam.width, ds.arr([3, 3, 3], "code_length"))
    assert_equal(cam.position, ds.arr([1, 1, 1], "code_length"))
    assert_equal(cam.focus, ds.arr([0, 0, 0], "code_length"))

    # test setting the attributes in various ways

    attribute_values = [
        (
            1,
            ds.arr([2, 2, 2], "code_length"),
        ),
        (
            [1],
            ds.arr([2, 2, 2], "code_length"),
        ),
        (
            [1, 2],
            RuntimeError,
        ),
        (
            [1, 1, 1],
            ds.arr([2, 2, 2], "code_length"),
        ),
        (
            (1, "code_length"),
            ds.arr([1, 1, 1], "code_length"),
        ),
        (
            ((1, "code_length"), (1, "code_length")),
            RuntimeError,
        ),
        (
            ((1, "cm"), (2, "cm"), (3, "cm")),
            ds.arr([0.5, 1, 1.5], "code_length"),
        ),
        (
            2 * u.cm,
            ds.arr([1, 1, 1], "code_length"),
        ),
        (
            ds.arr(2, "cm"),
            ds.arr([1, 1, 1], "code_length"),
        ),
        (
            [2 * u.cm],
            ds.arr([1, 1, 1], "code_length"),
        ),
        (
            [1, 2, 3] * u.cm,
            ds.arr([0.5, 1, 1.5], "code_length"),
        ),
        (
            [1, 2] * u.cm,
            RuntimeError,
        ),
        (
            [u.cm * w for w in [1, 2, 3]],
            ds.arr([0.5, 1, 1.5], "code_length"),
        ),
    ]

    # define default values to avoid accidentally setting focus = position
    default_values = {
        "focus": [0, 0, 0],
        "position": [4, 4, 4],
        "width": [1, 1, 1],
    }
    attribute_list = list(default_values.keys())

    for attribute in attribute_list:
        for other_attribute in [a for a in attribute_list if a != attribute]:
            setattr(cam, other_attribute, default_values[other_attribute])
        for attribute_value, expected_result in attribute_values:
            try:
                # test properties
                setattr(cam, attribute, attribute_value)
                assert_almost_equal(getattr(cam, attribute), expected_result)
            except RuntimeError:
                assert expected_result is RuntimeError

            try:
                # test setters/getters
                getattr(cam, f"set_{attribute}")(attribute_value)
                assert_almost_equal(getattr(cam, f"get_{attribute}")(), expected_result)
            except RuntimeError:
                assert expected_result is RuntimeError

    resolution_values = (
        (
            512,
            (512, 512),
        ),
        (
            (512, 512),
            (512, 512),
        ),
        (
            (256, 512),
            (256, 512),
        ),
        ((256, 256, 256), RuntimeError),
    )

    for resolution_value, expected_result in resolution_values:
        try:
            # test properties
            cam.resolution = resolution_value
            assert_equal(cam.resolution, expected_result)
        except RuntimeError:
            assert expected_result is RuntimeError

        try:
            # test setters/getters
            cam.set_resolution(resolution_value)
            assert_equal(cam.get_resolution(), expected_result)
        except RuntimeError:
            assert expected_result is RuntimeError

    for lens_type in valid_lens_types:
        cam.set_lens(lens_type)

    # See issue #1287
    cam.focus = [0, 0, 0]
    cam_pos = [1, 0, 0]
    north_vector = [0, 1, 0]
    cam.set_position(cam_pos, north_vector)
    cam_pos = [0, 1, 0]
    north_vector = [0, 0, 1]
    cam.set_position(cam_pos, north_vector)
