import numpy as np
from numpy.testing import assert_allclose, assert_raises

from yt.utilities.exceptions import YTException
from yt.utilities.orientation import (
    Orientation,
    _aligned,
    _validate_unit_vectors,
)
from yt.units.yt_array import YTArray


def test_aligned_detects_parallel_and_orthogonal_vectors():
    assert bool(_aligned(np.array([1.0, 0.0, 0.0]), np.array([-2.0, 0.0, 0.0])))
    assert not bool(
        _aligned(np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0]))
    )


def test_validate_unit_vectors_rejects_zero_and_aligned_vectors():
    normal_vector, north_vector = _validate_unit_vectors([0.0, 0.0, 2.0], [1.0, 0.0, 0.0])
    assert_allclose(normal_vector, [0.0, 0.0, 2.0])
    assert_allclose(north_vector, [1.0, 0.0, 0.0])

    with assert_raises(TypeError):
        _validate_unit_vectors(None, [1.0, 0.0, 0.0])

    with assert_raises(YTException):
        _validate_unit_vectors([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])

    with assert_raises(YTException):
        _validate_unit_vectors([0.0, 0.0, 1.0], [0.0, 0.0, 4.0])


def test_orientation_without_north_vector_builds_expected_basis():
    orientation = Orientation([0.0, 0.0, 1.0])

    assert orientation.steady_north is False
    assert_allclose(orientation.normal_vector, [0.0, 0.0, 1.0])
    assert_allclose(orientation.north_vector, [1.0, 0.0, 0.0])
    assert_allclose(
        orientation.unit_vectors,
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
    )
    assert_allclose(orientation.inv_mat, np.transpose(orientation.unit_vectors))


def test_orientation_for_x_axis_normal_builds_expected_axis_aligned_basis():
    orientation = Orientation([2.0, 0.0, 0.0])

    assert orientation.steady_north is False
    assert_allclose(orientation.normal_vector, [1.0, 0.0, 0.0])
    assert_allclose(orientation.north_vector, [0.0, 1.0, 0.0])
    assert_allclose(
        orientation.unit_vectors,
        [
            [0.0, 0.0, -1.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
        ],
    )


def test_orientation_with_explicit_north_vector_projects_and_normalizes():
    orientation = Orientation([0.0, 0.0, 2.0], [1.0, 1.0, 1.0])

    root_two = np.sqrt(2.0)
    assert orientation.steady_north is True
    assert_allclose(orientation.normal_vector, [0.0, 0.0, 1.0])
    assert_allclose(orientation.north_vector, [1.0 / root_two, 1.0 / root_two, 0.0])
    assert_allclose(
        orientation.unit_vectors,
        [
            [1.0 / root_two, -1.0 / root_two, 0.0],
            [1.0 / root_two, 1.0 / root_two, 0.0],
            [0.0, 0.0, 1.0],
        ],
    )
    assert_allclose(orientation.inv_mat, np.transpose(orientation.unit_vectors))


def test_setup_normalized_vectors_keeps_orthogonal_north_without_projection():
    orientation = Orientation.__new__(Orientation)
    orientation.steady_north = False

    orientation._setup_normalized_vectors(
        np.array([0.0, 0.0, 2.0]),
        np.array([0.0, 3.0, 0.0]),
    )

    assert_allclose(orientation.normal_vector, [0.0, 0.0, 1.0])
    assert_allclose(orientation.north_vector, [0.0, 1.0, 0.0])
    assert_allclose(orientation.unit_vectors[0], [1.0, 0.0, 0.0])


def test_orientation_init_uses_unit_vector_fallback_for_north_vector():
    class _FallbackOrientation(Orientation):
        def _setup_normalized_vectors(self, normal_vector, north_vector):
            self.unit_vectors = YTArray(
                [
                    [0.0, 1.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
                "",
            )
            self.north_vector = None

    orientation = _FallbackOrientation([0.0, 0.0, 1.0])

    assert_allclose(orientation.north_vector, [1.0, 0.0, 0.0])
