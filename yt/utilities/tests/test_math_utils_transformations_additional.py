import numpy as np
from numpy.testing import assert_allclose, assert_raises

from yt.units.yt_array import YTArray
from yt.utilities.math_utils import (
    compute_stddev_image,
    compute_cylindrical_radius,
    compute_parallel_velocity,
    compute_radial_velocity,
    compute_rotational_velocity,
    get_lookat_matrix,
    get_orthographic_matrix,
    get_perspective_matrix,
    get_rotation_matrix,
    get_scale_matrix,
    get_translate_matrix,
    modify_reference_frame,
    ortho_find,
    periodic_position,
    periodic_ray,
    quartiles,
    quaternion_mult,
    quaternion_to_rotation_matrix,
    rotate_vector_3D,
    rotation_matrix_to_quaternion,
)


class _FakePeriodicDataset:
    domain_left_edge = np.array([-1.0, -2.0, 0.0])
    domain_width = np.array([2.0, 4.0, 1.0])


def test_periodic_position_and_ray_have_exact_wrapped_outputs():
    wrapped = periodic_position(np.array([1.2, -5.5, 1.2]), _FakePeriodicDataset())
    assert_allclose(wrapped, [-0.8, -1.5, 0.2])

    segments = periodic_ray(
        YTArray([0.5, 0.5, 0.5], ""),
        YTArray([1.25, 1.25, 1.25], ""),
    )
    assert len(segments) == 2
    assert_allclose(segments[0][0], [0.5, 0.5, 0.5])
    assert_allclose(segments[0][1], [1.0, 1.0, 1.0])
    assert_allclose(segments[1][0], [0.0, 0.0, 0.0])
    assert_allclose(segments[1][1], [0.25, 0.25, 0.25])


def test_periodic_ray_with_explicit_bounds_handles_left_crossing():
    with np.errstate(invalid="ignore"):
        segments = periodic_ray(
            YTArray([0.1, 0.4, 0.5], ""),
            YTArray([-0.3, 0.4, 0.5], ""),
            left=YTArray([0.0, 0.0, 0.0], ""),
            right=YTArray([1.0, 1.0, 1.0], ""),
        )

    assert len(segments) == 2
    assert_allclose(segments[0][0], [0.1, 0.4, 0.5])
    assert_allclose(segments[0][1], [0.0, 0.4, 0.5])
    assert_allclose(segments[1][0], [1.0, 0.4, 0.5])
    assert_allclose(segments[1][1], [0.7, 0.4, 0.5])


def test_rotate_vector_3d_covers_axis_examples_and_input_validation():
    basis_vectors = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    assert_allclose(
        rotate_vector_3D(basis_vectors, 2, np.pi / 2.0),
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        atol=1.0e-15,
    )
    assert_allclose(
        rotate_vector_3D(np.array([1.0, 0.0, 0.0]), 1, np.pi / 2.0),
        [0.0, 0.0, 1.0],
        atol=1.0e-15,
    )
    assert_allclose(
        rotate_vector_3D(np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]), 0, np.pi / 2.0),
        [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0]],
        atol=1.0e-15,
    )

    with assert_raises(ValueError):
        rotate_vector_3D(np.zeros((2, 2)), 2, np.pi / 2.0)
    with assert_raises(ValueError):
        rotate_vector_3D(basis_vectors, 3, np.pi / 2.0)


def test_modify_reference_frame_handles_negative_z_and_general_rotation():
    negative_z_l = np.array([0.0, 0.0, -2.0])
    negative_z_p = np.array([[2.0, 1.0, 1.0], [1.0, 2.0, 1.0]])
    negative_z_v = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    l_out, p_out, v_out = modify_reference_frame(
        np.array([1.0, 1.0, 1.0]),
        negative_z_l,
        negative_z_p,
        negative_z_v,
    )
    assert_allclose(l_out, [0.0, 0.0, -2.0])
    assert_allclose(p_out, [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]])
    assert_allclose(v_out, [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]])

    l_rot, p_rot, v_rot = modify_reference_frame(
        np.array([0.5, 0.5, 0.5]),
        np.array([1.0, 0.0, 0.0]),
        np.array(
            [
                [1.0, 0.5, 0.5],
                [0.0, 0.5, 0.5],
                [0.5, 0.5, 0.5],
                [0.0, 0.0, 0.0],
            ]
        ),
        np.array(
            [
                [1.0, 0.5, 0.5],
                [0.0, 0.5, 0.5],
                [0.5, 0.5, 0.5],
                [0.0, 0.0, 0.0],
            ]
        ),
    )
    assert_allclose(l_rot, [0.0, 0.0, 1.0], atol=1.0e-14)
    assert_allclose(
        p_rot,
        [
            [0.0, 0.0, 0.5],
            [0.0, 0.0, -0.5],
            [0.0, 0.0, 0.0],
            [0.5, -0.5, -0.5],
        ],
        atol=1.0e-14,
    )
    assert_allclose(
        v_rot,
        [
            [-0.5, 0.5, 1.0],
            [-0.5, 0.5, 0.0],
            [-0.5, 0.5, 0.5],
            [0.0, 0.0, 0.0],
        ],
        atol=1.0e-14,
    )


def test_modify_reference_frame_position_only_and_velocity_only_paths():
    l_pos, p_only = modify_reference_frame(
        np.array([1.0, 1.0, 1.0]),
        np.array([0.0, 0.0, -2.0]),
        P=np.array([[2.0, 1.0, 1.0], [1.0, 2.0, 1.0]]),
        V=None,
    )
    assert_allclose(l_pos, [0.0, 0.0, -2.0])
    assert_allclose(p_only, [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]])

    l_aligned_vel, aligned_v_only = modify_reference_frame(
        np.array([1.0, 1.0, 1.0]),
        np.array([0.0, 0.0, -2.0]),
        P=None,
        V=np.array([[2.0, 1.0, 1.0], [1.0, 2.0, 1.0]]),
    )
    assert_allclose(l_aligned_vel, [0.0, 0.0, -2.0])
    assert_allclose(aligned_v_only, [[-2.0, -1.0, -1.0], [-1.0, -2.0, -1.0]])

    l_vel, v_only = modify_reference_frame(
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, -1.0, 0.0]),
        P=None,
        V=np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]),
    )
    assert_allclose(l_vel, [0.0, 0.0, np.sqrt(2.0)], atol=1.0e-14)
    assert_allclose(
        v_only,
        [
            [-1.0, 0.0, 0.0],
            [0.0, np.sqrt(0.5), np.sqrt(0.5)],
        ],
        atol=1.0e-14,
    )

    l_rot_pos, rotated_positions = modify_reference_frame(
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, -1.0, 0.0]),
        P=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        V=None,
    )
    assert_allclose(l_rot_pos, [0.0, 0.0, np.sqrt(2.0)], atol=1.0e-14)
    assert_allclose(
        rotated_positions,
        [
            [0.0, np.sqrt(0.5), np.sqrt(0.5)],
            [0.0, np.sqrt(0.5), -np.sqrt(0.5)],
        ],
        atol=1.0e-14,
    )


def test_velocity_and_radius_helpers_match_hand_computable_examples():
    center_of_mass = np.array([0.0, 0.0, 0.0])
    angular_momentum = np.array([0.0, 0.0, 1.0])
    positions = np.array(
        [
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
        ]
    )
    velocities = np.array(
        [
            [0.0, 1.0, 10.0],
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
        ]
    )

    with np.errstate(invalid="ignore"):
        rotational_velocity = compute_rotational_velocity(
            center_of_mass, angular_momentum, positions, velocities
        )
        radial_velocity = compute_radial_velocity(
            center_of_mass, angular_momentum, positions, velocities
        )

    assert_allclose(
        rotational_velocity,
        [1.0, 0.0, np.nan, np.sqrt(2.0)],
        atol=1.0e-15,
        equal_nan=True,
    )
    assert_allclose(
        compute_parallel_velocity(
            center_of_mass, angular_momentum, positions, velocities
        ),
        [10.0, -1.0, 1.0, -1.0],
    )
    assert_allclose(
        radial_velocity,
        [0.0, np.sqrt(2.0), np.nan, 0.0],
        atol=1.0e-15,
        equal_nan=True,
    )
    assert_allclose(
        compute_cylindrical_radius(
            center_of_mass, angular_momentum, positions, velocities
        ),
        [1.0, np.sqrt(2.0), 0.0, np.sqrt(2.0)],
        atol=1.0e-15,
    )


def test_ortho_find_covers_remaining_axis_selection_branches():
    with assert_raises(ValueError):
        ortho_find([0.0, 0.0, 0.0])

    vec1_y, vec2_y, vec3_y = ortho_find([1.0, 1.0, 0.0])
    assert_allclose(vec1_y, [np.sqrt(0.5), np.sqrt(0.5), 0.0], atol=1.0e-15)
    assert_allclose(vec2_y, [0.0, 0.0, 1.0], atol=1.0e-15)
    assert_allclose(vec3_y, [np.sqrt(0.5), -np.sqrt(0.5), 0.0], atol=1.0e-15)

    vec1_x, vec2_x, vec3_x = ortho_find([1.0, 0.0, 0.0])
    assert_allclose(vec1_x, [1.0, 0.0, 0.0], atol=1.0e-15)
    assert_allclose(vec2_x, [0.0, 1.0, 0.0], atol=1.0e-15)
    assert_allclose(vec3_x, [0.0, 0.0, 1.0], atol=1.0e-15)


def test_quartiles_covers_flattened_small_and_overwrite_examples():
    assert_allclose(quartiles(np.array([4.0, 1.0, 3.0, 2.0])), [1.5, 3.5])

    with assert_raises(AttributeError):
        quartiles(np.array([4.0, 1.0, 3.0, 2.0]), overwrite_input=True)

    two_by_two = np.array([[4.0, 1.0], [2.0, 3.0]])
    assert_allclose(
        quartiles(two_by_two.copy(), axis=0, overwrite_input=True),
        [[2.0, 1.0], [4.0, 3.0]],
    )

    rowwise = np.array([[3.0, 1.0, 2.0], [9.0, 7.0, 8.0]])
    out = np.empty((2,), dtype="float64")
    assert_allclose(quartiles(rowwise.copy(), axis=1, out=out), [[3.0, 9.0], [3.0, 9.0]])
    assert_allclose(
        quartiles(rowwise.copy(), axis=1, overwrite_input=True),
        [[1.0, 7.0], [3.0, 9.0]],
    )


def test_matrix_and_quaternion_helpers_return_expected_transforms():
    rotation_matrix = get_rotation_matrix(np.pi / 2.0, [0.0, 1.0, 0.0])
    assert_allclose(
        rotation_matrix,
        [
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
        ],
        atol=1.0e-15,
    )

    quaternion = rotation_matrix_to_quaternion(rotation_matrix)
    assert_allclose(quaternion, [np.sqrt(0.5), 0.0, -np.sqrt(0.5), 0.0], atol=1.0e-15)
    assert_allclose(quaternion_mult(np.array([1.0, 0.0, 0.0, 0.0]), quaternion), quaternion)
    assert_allclose(
        quaternion_to_rotation_matrix(quaternion),
        rotation_matrix,
        atol=1.0e-14,
    )

    assert_allclose(
        get_translate_matrix(1.0, 2.0, 3.0),
        [
            [1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 2.0],
            [0.0, 0.0, 1.0, 3.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
    )
    assert_allclose(
        get_scale_matrix(2.0, 3.0, 4.0),
        [
            [2.0, 0.0, 0.0, 0.0],
            [0.0, 3.0, 0.0, 0.0],
            [0.0, 0.0, 4.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
    )
    assert_allclose(
        get_lookat_matrix([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, -1.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
    )
    assert_allclose(
        get_perspective_matrix(90.0, 2.0, 1.0, 11.0),
        [
            [0.5, 0.0, 0.0, 0.0],
            [0.0, 2.0 / 3.0, -1.0 / 3.0, 0.0],
            [0.0, 0.0, -1.2, -2.2],
            [0.0, 0.0, -1.0, 0.0],
        ],
        atol=1.0e-7,
    )
    assert_allclose(
        get_orthographic_matrix(2.0, 2.0, 1.0, 11.0),
        [
            [0.25, 0.0, 0.0, 0.0],
            [0.0, 0.5, 0.0, 0.0],
            [0.0, 0.0, -0.2, 0.0],
            [0.0, 0.0, -1.2, 1.0],
        ],
        atol=1.0e-7,
    )


def test_rotation_matrix_to_quaternion_covers_x_y_and_z_dominant_branches():
    qx = rotation_matrix_to_quaternion(get_rotation_matrix(np.pi, [1.0, 0.0, 0.0]))
    qy = rotation_matrix_to_quaternion(get_rotation_matrix(np.pi, [0.0, 1.0, 0.0]))
    qz = rotation_matrix_to_quaternion(get_rotation_matrix(np.pi, [0.0, 0.0, 1.0]))

    assert_allclose(np.abs(qx), [0.0, 1.0, 0.0, 0.0], atol=1.0e-14)
    assert_allclose(np.abs(qy), [0.0, 0.0, 1.0, 0.0], atol=1.0e-14)
    assert_allclose(np.abs(qz), [0.0, 0.0, 0.0, 1.0], atol=1.0e-14)


def test_compute_stddev_image_handles_roundoff_and_large_negative_variances():
    assert_allclose(
        compute_stddev_image(np.array([1.0 - 1.0e-15, 4.0]), np.array([1.0, 2.0])),
        [0.0, 0.0],
        atol=1.0e-15,
    )

    with assert_raises(ValueError):
        compute_stddev_image(np.array([0.5, 4.0]), np.array([1.0, 2.0]))
