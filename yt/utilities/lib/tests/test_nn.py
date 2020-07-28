import numpy as np

from yt.testing import assert_array_equal
from yt.utilities.lib.bounded_priority_queue import (
    validate,
    validate_nblist,
    validate_pid,
)


# These test functions use utility functions in
# yt.utilities.lib.bounded_priority_queue
# to test functions which are not exposed at a python level
def test_bounded_priority_queue():
    dists = validate()
    answers = np.array([0.1, 0.001, -1.0, -1.0, -1.0])
    assert_array_equal(answers, dists)


def test_bounded_priority_queue_pid():
    dists, pids = validate_pid()
    answers = np.array([0.1, 0.001, -1.0, -1.0, -1.0])
    answers_pids = np.array([1, 10, -1, -1, -1])
    assert_array_equal(answers, dists)
    assert_array_equal(answers_pids, pids)


def test_neighbor_list():
    data, pids = validate_nblist()
    answers_data = np.array([1.0, 1.0, 1.0, 1.0])
    answers_pids = np.array([0, 1, 2, 3])
    assert_array_equal(answers_data, data)
    assert_array_equal(answers_pids, pids)
