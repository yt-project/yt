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
