import numpy as np
from numpy.testing import assert_allclose, assert_raises

from yt.utilities.math_utils import euclidean_dist, periodic_dist


def test_vectorized_periodic_distance_examples():
    a = np.array([[0.95, 0.25], [0.10, 0.40], [0.50, 0.90]])
    b = np.array([[0.05, 0.75], [0.10, 0.60], [0.50, 0.10]])

    assert_allclose(periodic_dist(a, b, 1.0), [0.1, 0.5744562646538028])
    assert_allclose(periodic_dist(b, a, 1.0), [0.1, 0.5744562646538028])
    assert_allclose(euclidean_dist(a, b), [0.9, 0.9643650760992956])
    assert_allclose(periodic_dist(a[:, 0], b[:, 0], 1.0), 0.1)


def test_mixed_periodicity_and_period_examples():
    a = [1.9, 0.1, 0.1]
    b = [0.1, 0.9, 0.5]

    assert_allclose(
        periodic_dist(a, b, [2.0, 1.0, 1.0], periodicity=(True, False, True)),
        0.9165151389911681,
    )
    assert_allclose(euclidean_dist(a, b), 2.009975124224178)


def test_periodic_distance_with_full_period_array_examples():
    a = np.array([[0.1, 0.4], [0.2, 0.3], [0.3, 0.1]])
    b = np.array([[0.9, 0.1], [0.2, 0.9], [0.1, 0.8]])
    period = np.array([[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]])

    assert_allclose(periodic_dist(a, b, period), [0.282842712474619, 0.58309518948453])


def test_periodic_distance_rejects_mismatched_shapes():
    with assert_raises(RuntimeError):
        periodic_dist([0.0, 0.0, 0.0], [[0.0], [0.0], [0.0]], 1.0)

    with assert_raises(ValueError):
        euclidean_dist([0.0, 0.0], [1.0, 1.0, 1.0])
