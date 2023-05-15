import re

import pytest
from unyt import unyt_array
from unyt.exceptions import UnitConversionError

from yt.testing import fake_amr_ds


@pytest.fixture(scope="module")
def reusable_fake_dataset():
    ds = fake_amr_ds(
        fields=[("gas", "density")],
        units=["g/cm**3"],
    )
    return ds


valid_single_str_values = ("center",)
valid_field_loc_str_values = ("min", "max")

DEFAUT_ERROR_MESSAGE = (
    "Expected any of the following\n"
    "- 'c', 'center', 'l', 'left', 'r', 'right', 'm', 'max', or 'min'\n"
    "- a 2 element tuple with 'min' or 'max' as the first element, followed by a field identifier\n"
    "- a 3 element array-like: for a unyt_array, expects length dimensions, otherwise code_lenght is assumed"
)


@pytest.mark.parametrize(
    "user_input",
    (
        # second element can be a single str or a field tuple (2 str), but not three
        (("max", ("not", "a", "field"))),
        # a 1-tuple is also not a valid field key
        (("max", ("notafield",))),
        # both elements need to be str
        (("max", (0, "invalid_field_type"))),
        (("max", ("invalid_field_type", 1))),
    ),
)
def test_invalid_center_type_default_error(reusable_fake_dataset, user_input):
    ds = reusable_fake_dataset
    with pytest.raises(
        TypeError,
        match=re.escape(f"Received {user_input!r}, ")
        + r"but failed to transform to a unyt_array \(obtained .+\)\.",
    ):
        # at the time of writing `axis` is an unused parameter of the base
        # sanitize center method, which is used directly for cartesian coordinate handlers
        # this probably hints that a refactor would make sense to separaet center sanitizing
        # and display_center calculation
        ds.coordinates.sanitize_center(user_input, axis=None)


@pytest.mark.parametrize(
    "user_input, error_type, error_message",
    (
        (
            "bad_str",
            ValueError,
            re.escape(
                "Received unknown center single string value 'bad_str'. "
                + DEFAUT_ERROR_MESSAGE
            ),
        ),
        (
            ("bad_str", ("gas", "density")),
            ValueError,
            re.escape(
                "Received unknown string value 'bad_str'. "
                f"Expected one of {valid_field_loc_str_values} (case insensitive)"
            ),
        ),
        (
            ("bad_str", "density"),
            ValueError,
            re.escape(
                "Received unknown string value 'bad_str'. "
                "Expected one of ('min', 'max') (case insensitive)"
            ),
        ),
        # even with exactly three elements, the dimension should be length
        (
            unyt_array([0.5] * 3, "kg"),
            UnitConversionError,
            "...",  # don't match the exact error message since it's unyt's responsibility
        ),
        # only validate 3 elements unyt_arrays
        (
            unyt_array([0.5] * 2, "cm"),
            TypeError,
            re.escape("Received unyt_array([0.5, 0.5], 'cm')"),
        ),
        (
            unyt_array([0.5] * 4, "cm"),
            TypeError,
            # don't attempt to match error message as details of how
            # a unyt array with more than a couple elements is displayed are out of our control
            "...",
        ),
        (
            # check that the whole shape is used in validation, not just the length (number of rows)
            unyt_array([0.5] * 6, "cm").reshape(3, 2),
            TypeError,
            # don't attempt to match error message as details of how
            # a unyt array with more than a couple elements is displayed are out of our control
            "...",
        ),
    ),
)
def test_invalid_center_special_cases(
    reusable_fake_dataset, user_input, error_type, error_message
):
    ds = reusable_fake_dataset
    with pytest.raises(error_type, match=error_message):
        ds.coordinates.sanitize_center(user_input, axis=None)
