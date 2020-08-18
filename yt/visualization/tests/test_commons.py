from pathlib import Path

from yt.visualization._commons import SUPPORTED_FORMATS, validate_image_name

filename_no_ext = Path("a_file_is_no_one")
filename_png = filename_no_ext.with_suffix(".png")
filename_pdf = filename_no_ext.with_suffix(".pdf")


def test_default():
    result = validate_image_name(filename_no_ext)
    expected = f"{filename_no_ext}.png"
    assert result == expected

    result = validate_image_name(filename_png)
    assert result == expected

    result = validate_image_name(filename_pdf)
    expected = str(filename_pdf)
    assert result == expected


def test_custom_valid_ext():
    for dext in SUPPORTED_FORMATS:
        ext = dext.replace(".", "")
        result1 = validate_image_name(filename_no_ext, suffix=ext)
        result2 = validate_image_name(filename_no_ext, suffix=dext)
        expected = f"{filename_no_ext}.{ext}"

        assert result1 == expected
        assert result2 == expected
