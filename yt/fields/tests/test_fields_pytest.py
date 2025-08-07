import pytest

from yt.fields.field_functions import validate_field_function


@pytest.mark.parametrize(
    "function",
    [
        pytest.param(lambda field, data: ..., id="field-data"),
        pytest.param(lambda data, field: ..., id="data-field"),
        pytest.param(lambda *, data, field: ..., id="data_kw-field_kw"),
        pytest.param(lambda data, *, field: ..., id="data-field_kw"),
        pytest.param(lambda field, *, data: ..., id="field-data_kw"),
        pytest.param(
            lambda data, field, extra=None: ..., id="data-field-extra_default"
        ),
    ],
)
def test_validate_field_function(function):
    assert validate_field_function(function) is None


@pytest.mark.parametrize(
    "function",
    [
        pytest.param(lambda data, field=None: ..., id="data-field_default"),
        pytest.param(lambda field, data=None: ..., id="field-data_default"),
    ],
)
def test_validate_field_function_default(function):
    with pytest.warns(
        UserWarning,
        match=r"These default values will never be used\. Drop them to avoid this warning\.$",
    ):
        validate_field_function(function)


@pytest.mark.parametrize(
    "function",
    [
        pytest.param(lambda: ..., id="no_args"),
        pytest.param(lambda dataa: ..., id="no-data"),
        pytest.param(lambda field, data, extra: ..., id="field-data-extra_nodefault"),
        pytest.param(lambda field, /, data: ..., id="field_pos-data"),
        pytest.param(lambda data, /, field: ..., id="data_pos-field"),
        pytest.param(lambda field, data, /: ..., id="field_pos-data_pos"),
        pytest.param(lambda *, data, field, extra: ..., id="data_kw-field_kw-extra_kw"),
        pytest.param(lambda data, field, *, extra: ..., id="data-field-extra_kw"),
        pytest.param(lambda data, *field: ..., id="data-field_varargs"),
        pytest.param(lambda field, *data: ..., id="field-data_varargs"),
        pytest.param(lambda data, **field: ..., id="data-field_kwargs"),
        pytest.param(lambda field, **data: ..., id="field-data_kwargs"),
    ],
)
def test_invalidate_field_function(function):
    with pytest.raises(
        TypeError,
        match=r"^Received field function .* with invalid signature",
    ):
        validate_field_function(function)
