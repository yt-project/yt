from functools import partial

from yt.data_objects.static_output import Dataset
from yt.testing import assert_raises
from yt.units import YTQuantity
from yt.units.unit_registry import UnitRegistry

mock_quan = partial(YTQuantity, registry=UnitRegistry())


def test_schema_validation():

    valid_schemas = [
        {"length_unit": 1.0},
        {"length_unit": [1.0]},
        {"length_unit": (1.0,)},
        {"length_unit": int(1.0)},
        {"length_unit": (1.0, "m")},
        {"length_unit": [1.0, "m"]},
        {"length_unit": YTQuantity(1.0, "m")},
    ]

    for schema in valid_schemas:
        uo = Dataset._sanitize_units_override(schema)
        for v in uo.values():
            q = mock_quan(v)  # check that no error (TypeError) is raised
            q.to("pc")  # check that q is a length


def test_invalid_schema_detection():
    invalid_key_schemas = [
        {"len_unit": 1.0},  # plain invalid key
        {"lenght_unit": 1.0},  # typo
    ]
    for invalid_schema in invalid_key_schemas:
        assert_raises(ValueError, Dataset._sanitize_units_override, invalid_schema)

    invalid_val_schemas = [
        {"length_unit": [1, 1, 1]},  # len(val) > 2
        {"length_unit": [1, 1, 1, 1, 1]},  # "data type not understood" in unyt
    ]

    for invalid_schema in invalid_val_schemas:
        assert_raises(TypeError, Dataset._sanitize_units_override, invalid_schema)

    # 0 shouldn't make sense
    invalid_number_schemas = [
        {"length_unit": 0},
        {"length_unit": [0]},
        {"length_unit": (0,)},
        {"length_unit": (0, "cm")},
    ]
    for invalid_schema in invalid_number_schemas:
        assert_raises(ValueError, Dataset._sanitize_units_override, invalid_schema)


def test_typing_error_detection():
    invalid_schema = {"length_unit": "1m"}

    # this is the error that is raised by unyt on bad input
    assert_raises(RuntimeError, mock_quan, invalid_schema["length_unit"])

    # check that the sanitizer function is able to catch the
    # type issue before passing down to unyt
    assert_raises(TypeError, Dataset._sanitize_units_override, invalid_schema)


def test_dimensionality_error_detection():
    invalid_schema = {"length_unit": YTQuantity(1.0, "s")}
    assert_raises(ValueError, Dataset._sanitize_units_override, invalid_schema)
