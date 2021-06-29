from nose.tools import assert_raises

from yt.funcs import just_one, levenshtein_distance, validate_axis, validate_center
from yt.testing import assert_equal, fake_amr_ds
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
    desired = (
        "Expected axis of int or char type (can be "
        "[0, 'x', 'X', 1, 'y', 'Y', 2, 'z', 'Z']), received 'r'."
    )
    assert_equal(str(ex.exception)[:40], desired[:40])


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


def test_levenshtein():
    assert_equal(levenshtein_distance("abcdef", "abcdef"), 0)

    # Deletions / additions
    assert_equal(levenshtein_distance("abcdef", "abcde"), 1)
    assert_equal(levenshtein_distance("abcdef", "abcd"), 2)
    assert_equal(levenshtein_distance("abcdef", "abc"), 3)

    assert_equal(levenshtein_distance("abcdf", "abcdef"), 1)
    assert_equal(levenshtein_distance("cdef", "abcdef"), 2)
    assert_equal(levenshtein_distance("bde", "abcdef"), 3)

    # Substitutions
    assert_equal(levenshtein_distance("abcd", "abc_"), 1)
    assert_equal(levenshtein_distance("abcd", "ab__"), 2)
    assert_equal(levenshtein_distance("abcd", "a___"), 3)
    assert_equal(levenshtein_distance("abcd", "____"), 4)

    # Deletion + Substitutions
    assert_equal(levenshtein_distance("abcd", "abc_z"), 2)
    assert_equal(levenshtein_distance("abcd", "ab__zz"), 4)
    assert_equal(levenshtein_distance("abcd", "a___zzz"), 6)
    assert_equal(levenshtein_distance("abcd", "____zzzz"), 8)

    # Max distance
    assert_equal(levenshtein_distance("abcd", "", max_dist=0), 1)
    assert_equal(levenshtein_distance("abcd", "", max_dist=3), 4)
    assert_equal(levenshtein_distance("abcd", "", max_dist=10), 4)
    assert_equal(levenshtein_distance("abcd", "", max_dist=1), 2)
    assert_equal(levenshtein_distance("abcd", "a", max_dist=2), 3)
    assert_equal(levenshtein_distance("abcd", "ad", max_dist=2), 2)
    assert_equal(levenshtein_distance("abcd", "abd", max_dist=2), 1)
    assert_equal(levenshtein_distance("abcd", "abcd", max_dist=2), 0)
