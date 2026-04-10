import numpy as np
from numpy.testing import assert_array_equal

from yt.utilities.nodal_data_utils import (
    _get_indices,
    _get_linear_index,
    get_nodal_data,
    get_nodal_slices,
)


class _DummyFieldInfo:
    def __init__(self, nodal_flag):
        self.nodal_flag = nodal_flag


class _DummyDS:
    def __init__(self, nodal_flag):
        self.nodal_flag = nodal_flag

    def _get_field_info(self, field):
        return _DummyFieldInfo(self.nodal_flag)


class _DummyDataSource(dict):
    def __init__(self, nodal_flag, field, field_data):
        super().__init__({field: field_data})
        self.ds = _DummyDS(nodal_flag)


def test_get_linear_index_and_indices_lookup():
    assert _get_linear_index(np.array([0, 0])) == 0
    assert _get_linear_index(np.array([0, 1])) == 1
    assert _get_linear_index(np.array([1, 0])) == 2
    assert _get_linear_index(np.array([1, 1])) == 3

    assert _get_linear_index(np.array([1, 0, 1])) == 5
    assert _get_linear_index(np.array([1, 1, 1])) == 7

    assert_array_equal(_get_indices(np.array([1, 0])), [0, 0, 1, 1, 0, 0, 1, 1])
    assert_array_equal(_get_indices(np.array([1, 0, 1])), [0, 1, 0, 1, 2, 3, 2, 3])


def test_get_nodal_data_selects_expected_columns():
    field = ("gas", "test")
    field_data = np.arange(16).reshape(2, 8)

    data_source_2d = _DummyDataSource(np.array([1, 0]), field, field_data)
    expected_2d = field_data[:, [0, 0, 1, 1, 0, 0, 1, 1]]
    assert_array_equal(get_nodal_data(data_source_2d, field), expected_2d)

    data_source_3d = _DummyDataSource(np.array([1, 1, 1]), field, field_data)
    assert_array_equal(get_nodal_data(data_source_3d, field), field_data)


def test_get_nodal_slices_in_2d_and_3d():
    slices_2d = get_nodal_slices((4, 5), np.array([0, 1]), 2)
    expected_2d = (
        [slice(0, 4), slice(0, 4)],
        [slice(0, 4), slice(1, 5)],
    )
    assert slices_2d == expected_2d

    slices_3d = get_nodal_slices((4, 5, 6), np.array([1, 0, 1]), 3)
    expected_3d = (
        [slice(0, 3), slice(0, 5), slice(0, 5)],
        [slice(0, 3), slice(0, 5), slice(1, 6)],
        [slice(1, 4), slice(0, 5), slice(0, 5)],
        [slice(1, 4), slice(0, 5), slice(1, 6)],
    )
    assert slices_3d == expected_3d
