from functools import partial

import pytest

from yt import derived_field
from yt.fields import local_fields
from yt.testing import fake_random_ds


def test_add_field_lambda():
    ds = fake_random_ds(16)

    ds.add_field(
        ("gas", "spam"),
        lambda field, data: data["gas", "density"],
        sampling_type="cell",
    )

    # check access
    ds.all_data()["gas", "spam"]


def test_add_field_partial():
    ds = fake_random_ds(16)

    def _spam(field, data, factor):
        return factor * data["gas", "density"]

    ds.add_field(
        ("gas", "spam"),
        partial(_spam, factor=1),
        sampling_type="cell",
    )

    # check access
    ds.all_data()["gas", "spam"]


def test_add_field_arbitrary_callable():
    ds = fake_random_ds(16)

    class Spam:
        def __call__(self, field, data):
            return data["gas", "density"]

    ds.add_field(("gas", "spam"), Spam(), sampling_type="cell")

    # check access
    ds.all_data()["gas", "spam"]


def test_add_field_uncallable():
    ds = fake_random_ds(16)

    class Spam:
        pass

    with pytest.raises(TypeError, match=r"(is not a callable object)$"):
        ds.add_field(("bacon", "spam"), Spam(), sampling_type="cell")


def test_add_field_wrong_signature():
    ds = fake_random_ds(16)

    def _spam(data, field):
        return data["gas", "density"]

    with pytest.raises(
        TypeError,
        match=(
            r"Received field function <function .*> with invalid signature\. "
            r"Expected exactly 2 positional parameters \('field', 'data'\), got \('data', 'field'\)"
        ),
    ):
        ds.add_field(("bacon", "spam"), _spam, sampling_type="cell")


def test_add_field_keyword_only():
    ds = fake_random_ds(16)

    def _spam(field, *, data):
        return data["gas", "density"]

    with pytest.raises(
        TypeError,
        match=(
            r"Received field function .* with invalid signature\. "
            r"Parameters 'field' and 'data' must accept positional values \(they cannot be keyword-only\)"
        ),
    ):
        ds.add_field(
            ("bacon", "spam"),
            _spam,
            sampling_type="cell",
        )


def test_derived_field(monkeypatch):

    tmp_field_info = local_fields.LocalFieldInfoContainer(None, [], None)
    monkeypatch.setattr(local_fields, "local_fields", tmp_field_info)

    @derived_field(name="pressure", sampling_type="cell", units="dyne/cm**2")
    def _pressure(field, data):
        return (
            (data.ds.gamma - 1.0)
            * data["gas", "density"]
            * data["gas", "specific_thermal_energy"]
        )
