import os

from nose.tools import assert_raises
from numpy.testing import assert_equal

from yt.funcs import just_one, simple_download_file, validate_axis, validate_center
from yt.testing import fake_amr_ds
from yt.units import YTArray, YTQuantity


def test_validate_axis():
    validate_axis(None, 0)
    validate_axis(None, "X")

    ds = fake_amr_ds(geometry="cylindrical")
    ds.slice("Theta", 0.25)

    with assert_raises(TypeError) as ex:
        # default geometry is cartesian
        ds = fake_amr_ds()
        ds.slice("r", 0.25)
    desired = "Expected axis to be any of [0, 1, 2, 'x', 'y', 'z', 'X', 'Y', 'Z'], received 'r'"

    actual = str(ex.exception)
    assert actual == desired


def test_validate_center():
    validate_center("max")
    validate_center("MIN_")

    with assert_raises(TypeError) as ex:
        validate_center("avg")
    desired = (
        "Expected 'center' to be in ['c', 'center', 'm', 'max', 'min'] "
        "or the prefix to be 'max_'/'min_', received 'avg'."
    )
    assert_equal(str(ex.exception), desired)

    validate_center(YTQuantity(0.25, "cm"))
    validate_center([0.25, 0.25, 0.25])

    class CustomCenter:
        def __init__(self, center):
            self.center = center

    with assert_raises(TypeError) as ex:
        validate_center(CustomCenter(10))
    desired = (
        "Expected 'center' to be a numeric object of type "
        "list/tuple/np.ndarray/YTArray/YTQuantity, received "
        "'yt.tests.test_funcs.test_validate_center.<locals>."
        "CustomCenter'."
    )
    assert_equal(str(ex.exception)[:50], desired[:50])


def test_just_one():
    # Check that behaviour of this function is consistent before and after refactor
    # PR 2893
    for unit in ["mm", "cm", "km", "pc", "g", "kg", "M_sun"]:
        obj = YTArray([0.0, 1.0], unit)
        expected = YTQuantity(obj.flat[0], obj.units, registry=obj.units.registry)
        jo = just_one(obj)
        assert jo == expected


def test_simple_download_file():
    fn = simple_download_file("http://yt-project.org", "simple-download-file")
    try:
        assert fn == "simple-download-file"
        assert os.path.exists("simple-download-file")
    finally:
        # Clean up after ourselves.
        try:
            os.unlink("simple-download-file")
        except FileNotFoundError:
            pass

    with assert_raises(RuntimeError) as ex:
        simple_download_file("http://yt-project.org/404", "simple-download-file")

    desired = "Attempt to download file from http://yt-project.org/404 failed with error 404: Not Found."
    actual = str(ex.exception)
    assert actual == desired
    assert not os.path.exists("simple-download-file")
