import pytest
from numpy.testing import assert_raises

from yt.visualization._commons import (
    _swap_arg_pair_order,
    _swap_axes_extents,
    validate_image_name,
)


@pytest.mark.parametrize(
    "name, expected",
    [
        ("noext", "noext.png"),
        ("nothing.png", "nothing.png"),
        ("nothing.pdf", "nothing.pdf"),
        ("version.1.2.3", "version.1.2.3.png"),
    ],
)
def test_default(name, expected):
    result = validate_image_name(name)
    assert result == expected


@pytest.mark.parametrize(
    "name, suffix, expected",
    [
        ("noext", ".png", "noext.png"),
        ("noext", None, "noext.png"),
        ("nothing.png", ".png", "nothing.png"),
        ("nothing.png", None, "nothing.png"),
        ("nothing.png", ".pdf", "nothing.pdf"),
        ("nothing.pdf", ".pdf", "nothing.pdf"),
        ("nothing.pdf", None, "nothing.pdf"),
        ("nothing.pdf", ".png", "nothing.png"),
        ("version.1.2.3", ".png", "version.1.2.3.png"),
        ("version.1.2.3", None, "version.1.2.3.png"),
        ("version.1.2.3", ".pdf", "version.1.2.3.pdf"),
    ],
)
@pytest.mark.filterwarnings(
    r"ignore:Received two valid image formats '\w+' \(from filename\) "
    r"and '\w+' \(from suffix\). The former is ignored.:UserWarning"
)
def test_custom_valid_ext(name, suffix, expected):
    result1 = validate_image_name(name, suffix=suffix)
    assert result1 == expected
    if suffix is not None:
        alt_suffix = suffix.replace(".", "")
        result2 = validate_image_name(name, suffix=alt_suffix)
        assert result2 == expected


def test_extent_swap():
    input_extent = [1, 2, 3, 4]
    expected = [3, 4, 1, 2]
    assert _swap_axes_extents(input_extent) == expected
    assert _swap_axes_extents(tuple(input_extent)) == tuple(expected)


def test_swap_arg_pair_order():
    assert _swap_arg_pair_order(1, 2) == (2, 1)
    assert _swap_arg_pair_order(1, 2, 3, 4, 5, 6) == (2, 1, 4, 3, 6, 5)
    assert_raises(TypeError, _swap_arg_pair_order, 1)
    assert_raises(TypeError, _swap_arg_pair_order, 1, 2, 3)
