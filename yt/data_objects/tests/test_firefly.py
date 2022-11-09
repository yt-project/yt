import tempfile
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from yt.testing import fake_particle_ds, requires_module
from yt.utilities.exceptions import YTFieldNotFound


@requires_module("Firefly")
def test_firefly_JSON_string():

    ds = fake_particle_ds()
    ad = ds.all_data()
    reader = ad.create_firefly_object(
        None,
        velocity_units="cm/s",
        coordinate_units="cm",
    )

    reader.dumpToJSON(write_jsons_to_disk=False)

    ## reader.JSON was not output to string correctly
    ##  either Firefly is damaged or needs a hotfix-- try reinstalling.
    ##  if that doesn't work contact the developers
    ##  at github.com/ageller/Firefly/issues.
    assert len(reader.JSON) > 0


@requires_module("Firefly")
def test_firefly_write_to_disk():
    tmpdir = tempfile.mkdtemp()

    ds = fake_particle_ds()
    ad = ds.all_data()
    reader = ad.create_firefly_object(
        tmpdir,
        velocity_units="cm/s",
        coordinate_units="cm",
    )

    reader.dumpToJSON()


@pytest.fixture
def firefly_test_dataset():
    # create dataset
    ds_fields = [
        # Assumed present
        ("PartType1", "particle_position_x"),
        ("PartType1", "particle_position_y"),
        ("PartType1", "particle_position_z"),
        ("gas", "particle_position_x"),
        ("gas", "particle_position_y"),
        ("gas", "particle_position_z"),
        # User input
        ("PartType1", "Masses"),
        ("gas", "Masses"),
        ("gas", "Temperature"),
    ]
    ds_field_units = [
        "code_length",
        "code_length",
        "code_length",
        "code_length",
        "code_length",
        "code_length",
        "code_mass",
        "code_mass",
        "code_temperature",
    ]
    ds_negative = [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]
    ds = fake_particle_ds(
        fields=ds_fields,
        units=ds_field_units,
        negative=ds_negative,
    )
    dd = ds.all_data()
    return (ds, dd)


@requires_module("firefly")
@pytest.mark.parametrize(
    "fields_to_include,fields_units,pgs_to_test,expectation",
    [
        (None, None, None, does_not_raise()),  # Test default values
        ([], [], None, does_not_raise()),  # Test empty fields
        (
            ["Masses"],
            ["code_mass"],
            None,
            does_not_raise(),
        ),  # Test common field (Masses)
        (
            ["Temperature"],
            ["code_temperature"],
            "gas",
            does_not_raise(),
        ),  # Test unique field (Temperature)
        (
            ["dinos"],
            ["code_length"],
            None,
            pytest.raises(YTFieldNotFound),
        ),  # Test nonexistent field (dinos)
        (
            [("gas", "Temperature")],
            ["code_temperature"],
            None,
            does_not_raise(),
        ),  # Test existing field tuple (gas, Temperature)
        (
            [("PartType1", "Temperature")],
            ["code_temperature"],
            None,
            pytest.raises(YTFieldNotFound),
        ),  # Test nonexistent field tuple (PartType1, Temperature)
    ],
)
def test_fields_specification(
    firefly_test_dataset, fields_to_include, fields_units, pgs_to_test, expectation
):
    if pgs_to_test is None:
        pgs_to_test = ["PartType1", "gas"]

    (ds, dd) = firefly_test_dataset
    with expectation:
        reader = dd.create_firefly_object(
            fields_to_include=fields_to_include,
            fields_units=fields_units,
            coordinate_units="code_length",
        )
    if not isinstance(expectation, does_not_raise):
        # constructed like this so we don't need to import RaiseError from pytest
        return
    assert_array_equal(
        dd[("PartType1", "relative_particle_position")].d,
        reader.particleGroups[0].coordinates,
    )
    assert_array_equal(
        dd[("gas", "relative_particle_position")].d,
        reader.particleGroups[1].coordinates,
    )
    if fields_to_include is None:
        return
    all_pgs = ["PartType1", "gas"]
    for field in fields_to_include:
        if isinstance(field, tuple):
            pg = reader.particleGroups[0 if field[0] == "PartType1" else 1]
            assert field[1] in pg.field_names
            arrind = np.flatnonzero(pg.field_names == field[1])[0]
            assert_array_equal(dd[field].d, pg.field_arrays[arrind])
        else:
            for ptype in all_pgs:
                pg = reader.particleGroups[0 if ptype == "PartType1" else 1]
                if ptype not in pgs_to_test:
                    assert field not in pg.field_names
                    continue
                assert field in pg.field_names
                arrind = np.flatnonzero(pg.field_names == field)[0]
                assert_array_equal(dd[ptype, field].d, pg.field_arrays[arrind])
