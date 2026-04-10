import numpy as np
from numpy.testing import assert_array_almost_equal

from yt.utilities.math_utils import (
    get_cyl_r,
    get_cyl_r_component,
    get_cyl_theta,
    get_cyl_theta_component,
    get_cyl_z,
    get_cyl_z_component,
    get_sph_phi,
    get_sph_phi_component,
    get_sph_r,
    get_sph_r_component,
    get_sph_theta,
    get_sph_theta_component,
)


def test_axis_aligned_coordinate_conversion_examples():
    coords = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ]
    )
    normal = [0.0, 0.0, 1.0]

    assert_array_almost_equal(get_sph_r(coords), [1.0, 1.0, 1.0, 0.0])
    assert_array_almost_equal(get_sph_theta(coords, normal), [np.pi / 2, np.pi / 2, 0.0, 0.0])
    assert_array_almost_equal(get_sph_phi(coords, normal), [0.0, np.pi / 2, 0.0, 0.0])

    assert_array_almost_equal(get_cyl_r(coords, normal), [1.0, 1.0, 0.0, 0.0])
    assert_array_almost_equal(get_cyl_theta(coords, normal), [0.0, np.pi / 2, 0.0, 0.0])
    assert_array_almost_equal(get_cyl_z(coords, normal), [0.0, 0.0, 1.0, 0.0])


def test_axis_aligned_component_examples():
    normal = [0.0, 0.0, 1.0]
    theta = np.array([np.pi / 2])
    phi = np.array([0.0])
    vectors = np.array([[3.0], [4.0], [5.0]])

    assert_array_almost_equal(get_cyl_r_component(vectors, phi, normal), [3.0])
    assert_array_almost_equal(get_cyl_theta_component(vectors, phi, normal), [4.0])
    assert_array_almost_equal(get_cyl_z_component(vectors, normal), [5.0])

    assert_array_almost_equal(get_sph_r_component(vectors, theta, phi, normal), [3.0])
    assert_array_almost_equal(get_sph_phi_component(vectors, phi, normal), [4.0])
    assert_array_almost_equal(get_sph_theta_component(vectors, theta, phi, normal), [-5.0])
