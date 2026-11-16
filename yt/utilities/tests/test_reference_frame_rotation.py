import numpy as np
import pytest
from numpy.testing import assert_allclose

from yt.utilities.math_utils import modify_reference_frame


@pytest.mark.parametrize("omega_z", [-1.0, -0.25, -7.5])
@pytest.mark.parametrize("n_particles", [1, 5, 100])
def test_modify_reference_frame_L_along_negative_z(n_particles, omega_z):
    """
    Checks that L is rotated onto +z.

    When L lies exactly along -z, modify_reference_frame
    takes a special case where L is rotated 180 degrees
    about the x-axis so that it points along +z.

    Particles are generated at random positions in the xy-plane
    and given velocities for a clockwise rotation about z. By keeping
    z = 0, this makes sure that the total L lies along -z so that the
    special case is triggered.

    The returned L must point along +z with the same magnitude as the original L.
    Lastly, the angular momentum recomputed from the returned P and V must match L.
    """
    rng = np.random.default_rng(seed=0)
    P = rng.uniform(-10, 10, size=(n_particles, 3))
    P[:, 2] = 0.0

    omega = np.array([0.0, 0.0, omega_z])

    V = np.cross(omega, P)
    L = np.sum(np.cross(P, V), axis=0)

    L_expected = np.array([0.0, 0.0, np.linalg.norm(L)])

    L_new, P_new, V_new = modify_reference_frame(np.zeros(3), L, P, V)

    atol = 1e-12 * np.linalg.norm(L)

    assert isinstance(L_new, np.ndarray)
    assert_allclose(L_new, L_expected, atol=atol)
    assert_allclose(np.sum(np.cross(P_new, V_new), axis=0), L_expected, atol=atol)
