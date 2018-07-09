import numpy as np

from yt.utilities.lib.bounded_priority_queue import validate, validate_pid

from yt.testing import \
    fake_random_ds, \
    assert_array_equal, \
    assert_array_less, \
    assert_equal, \
    assert_raises

def test_bounded_priority_queue():
    dists = validate()
    answers = np.array([0.1, 0.001, -1., -1., -1.])
    assert_array_equal(answers, dists)

def test_bounded_priority_queue_pid():
    dists, pids = validate_pid()
    answers = np.array([0.1, 0.001, -1., -1., -1.])
    answers_pids = np.array([ 1, 10, -1, -1, -1])
    assert_array_equal(answers, dists)
    assert_array_equal(answers_pids, pids)

