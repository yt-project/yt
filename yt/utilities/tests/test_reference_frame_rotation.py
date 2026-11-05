import numpy as np
from numpy.testing import assert_allclose

from yt.utilities.math_utils import modify_reference_frame


def test_modify_reference_frame():
    # Three particles circling clockwise about +z
    P = np.array([[2.0, 0.0, 0.0], [0.0, -3.0, 0.0], [-1.0, 0.0, 0.0]])
    V = np.array([[0.0, -1.0, 0.0], [-2.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
    L = np.sum(np.cross(P, V), axis=0)

    L_new, P_new, V_new = modify_reference_frame(np.zeros(3), L, P.copy(), V.copy())

    assert_allclose(L_new, [0.0, 0.0, 12.0], atol=1e-12)
    assert_allclose(
        np.sum(np.cross(P_new, V_new), axis=0), [0.0, 0.0, 12.0], atol=1e-12
    )
