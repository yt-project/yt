import logging
from functools import partial

import numpy as np
import pytest
import unyt

import yt
from yt import derived_field
from yt.fields import local_fields
from yt.testing import fake_random_ds
from yt.utilities.logger import ytLogger


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


def test_add_field_flipped_signature():
    # before yt 4.5, the only valid signature was
    # `function(data)`, but now we allow `function(data, field)`
    ds = fake_random_ds(16)

    def _spam(data, field):
        return data["gas", "density"]

    ds.add_field(("bacon", "spam"), _spam, sampling_type="cell")


def test_add_field_signature_v2():
    # before yt 4.5, the only valid signature was
    # `function(data)`, but now we allow `function(data)`
    ds = fake_random_ds(16)

    def _spam(data):
        return data["gas", "density"]

    ds.add_field(("bacon", "spam"), _spam, sampling_type="cell")


def test_add_field_keyword_only():
    ds = fake_random_ds(16)

    def _spam(field, *, data):
        return data["gas", "density"]

    ds.add_field(
        ("bacon", "spam"),
        _spam,
        sampling_type="cell",
    )


def test_derived_field(monkeypatch):
    tmp_field_info = local_fields.LocalFieldInfoContainer(None, [], None)
    monkeypatch.setattr(local_fields, "local_fields", tmp_field_info)

    @derived_field(name="pressure", sampling_type="cell", units="dyne/cm**2")
    def _pressure(data):
        return (
            (data.ds.gamma - 1.0)
            * data["gas", "density"]
            * data["gas", "specific_thermal_energy"]
        )


@pytest.mark.parametrize(
    "add_field_kwargs",
    [
        # full default: auto unit detection, no (in)validation
        {},
        # explicit "auto", should be identical to default behaviour
        {"units": "auto"},
        # explicitly requesting dimensionless units
        {"units": "dimensionless"},
        # explicitly requesting dimensionless units (short hand)
        {"units": ""},
        # explictly requesting no dimensions
        {"dimensions": yt.units.dimensionless},
        # should work with unyt.dimensionless too
        {"dimensions": unyt.dimensionless},
        # supported short hand
        {"dimensions": "dimensionless"},
    ],
)
def test_dimensionless_field(add_field_kwargs):
    ds = fake_random_ds(16)

    def _dimensionless_field(data):
        return data["gas", "density"] / data["gas", "density"].units

    ds.add_field(
        name=("gas", "dimensionless_density"),
        function=_dimensionless_field,
        sampling_type="local",
        **add_field_kwargs,
    )
    # check access
    ds.all_data()["gas", "dimensionless_density"]


def test_add_field_quick():
    ds = fake_random_ds(16)

    # Test subtractions
    for field in (
        -ds.fields.gas.density,
        -1 * ds.fields.gas.density,
        ds.fields.gas.density - 2 * ds.fields.gas.density,
        ds.fields.gas.density * (-1),
    ):
        ds.add_field(("gas", "-density"), field, force_override=True)
        np.testing.assert_allclose(ds.r["gas", "-density"], -ds.r["gas", "density"])

    # Test additions
    for field in (
        ds.fields.index.ones + 1,
        1 + ds.fields.index.ones,
    ):
        ds.add_field(("index", "two"), field, force_override=True)
        np.testing.assert_allclose(ds.r["index", "two"], 2)

    # Test multiplications
    for field in (ds.fields.gas.density * 2, 2 * ds.fields.gas.density):
        ds.add_field(("gas", "twodensity"), field, force_override=True)
        np.testing.assert_allclose(
            ds.r["gas", "twodensity"], ds.r["gas", "density"] * 2
        )

    # Test divisions
    for field in (ds.fields.gas.density / 2, 0.5 * ds.fields.gas.density):
        ds.add_field(("gas", "halfdensity"), field, force_override=True)
        np.testing.assert_allclose(
            ds.r["gas", "halfdensity"], ds.r["gas", "density"] / 2
        )

    # Test right division
    for field in (1 / ds.fields.gas.density,):
        ds.add_field(("gas", "one_over_density"), field, force_override=True)
        np.testing.assert_allclose(
            ds.r["gas", "one_over_density"], 1 / ds.r["gas", "density"]
        )


def test_add_field_quick_syntax2():
    fields = ("density", "temperature")
    units = ("g/cm**3", "K")
    ds = fake_random_ds(16, fields=fields, units=units)

    # Left multiplication
    ds.fields.gas.density_two = ds.fields.gas.density * 2
    np.testing.assert_allclose(ds.r["gas", "density_two"], ds.r["gas", "density"] * 2)
    ds.fields.gas.density_two2 = ds.fields.gas.density + ds.fields.gas.density
    np.testing.assert_allclose(ds.r["gas", "density_two2"], ds.r["gas", "density"] * 2)

    # Right multiplication
    ds.fields.gas.two_density = 2 * ds.fields.gas.density
    np.testing.assert_allclose(ds.r["gas", "two_density"], ds.r["gas", "density"] * 2)

    # Left division
    ds.fields.gas.half_density = ds.fields.gas.density / 2
    np.testing.assert_allclose(ds.r["gas", "half_density"], ds.r["gas", "density"] / 2)
    ds.fields.gas.density_over_temperature = (
        ds.fields.gas.density / ds.fields.gas.temperature
    )
    np.testing.assert_allclose(
        ds.r["gas", "density_over_temperature"],
        ds.r["gas", "density"] / ds.r["gas", "temperature"],
    )

    # Right division
    ds.fields.gas.one_over_density = 1 / ds.fields.gas.density
    np.testing.assert_allclose(
        ds.r["gas", "one_over_density"], 1 / ds.r["gas", "density"]
    )

    # Subtraction
    ds.fields.gas.neg_density = -ds.fields.gas.density
    np.testing.assert_allclose(ds.r["gas", "neg_density"], -ds.r["gas", "density"])
    ds.fields.gas.density_minus_twice_density = (
        ds.fields.gas.density - 2 * ds.fields.gas.density
    )
    np.testing.assert_allclose(
        ds.r["gas", "density_minus_twice_density"],
        -ds.r["gas", "density"],
    )

    # Complex expression
    ds.fields.gas.kbT_per_V = (
        ds.fields.gas.temperature * ds.units.kb / ds.fields.gas.volume
    )
    np.testing.assert_allclose(
        ds.r["gas", "kbT_per_V"],
        ds.r["gas", "temperature"] * ds.units.kb / ds.r["gas", "volume"],
    )

    # Returning boolean
    dx_min = ds.r["index", "dx"].min()
    ds.fields.gas.smallest_cells = ds.fields.gas.dx == dx_min
    np.testing.assert_allclose(
        ds.r["gas", "smallest_cells"].value, (dx_min == ds.r["gas", "dx"])
    )


@pytest.fixture()
def capturable_logger(caplog):
    """
    This set the minimal conditions to make pytest's caplog fixture usable.
    """

    propagate = ytLogger.propagate
    ytLogger.propagate = True

    with caplog.at_level(logging.WARNING, "yt"):
        yield

    ytLogger.propagate = propagate


@pytest.mark.usefixtures("capturable_logger")
def test_add_field_quick_syntax_warnings(caplog):
    # Make sure we get a warning when overriding
    # an existing field with a different function
    ds = fake_random_ds(16)

    warn_str = "Field ('gas', 'density2') already exists. To override use `force_override=True`"
    # First time, no warning
    ds.add_field(("gas", "density2"), ds.fields.gas.density * 2)
    assert warn_str not in caplog.text

    # Second time, warning
    caplog.clear()
    ds.add_field(("gas", "density2"), ds.fields.gas.density * 2)
    assert warn_str in caplog.text

    warn_str = "Field ('gas', 'density3') already exists. To override use `force_override=True`"
    # Third time, new definition, no warning
    caplog.clear()
    ds.fields.gas.density3 = ds.fields.gas.density * 3
    assert warn_str not in caplog.text
    caplog.clear()

    # Fourth time: new definition, warning!
    caplog.clear()
    ds.fields.gas.density3 = ds.fields.gas.density * 3
    assert warn_str in caplog.text
