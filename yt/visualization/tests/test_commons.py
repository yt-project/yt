import pytest

from yt.visualization._commons import validate_image_name


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
        ("nothing.png", ".png", "nothing.png"),
        ("nothing.png", ".pdf", "nothing.pdf"),
        ("nothing.pdf", ".pdf", "nothing.pdf"),
        ("nothing.pdf", ".png", "nothing.png"),
        ("version.1.2.3", ".png", "version.1.2.3.png"),
        ("version.1.2.3", ".pdf", "version.1.2.3.pdf"),
    ],
)
def test_custom_valid_ext(name, suffix, expected):
    alt_suffix = suffix.replace(".", "")
    result1 = validate_image_name(name, suffix=suffix)
    result2 = validate_image_name(name, suffix=alt_suffix)

    assert result1 == expected
    assert result2 == expected
