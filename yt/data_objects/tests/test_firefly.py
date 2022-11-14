import tempfile

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from yt.testing import fake_particle_ds, requires_module
from yt.utilities.exceptions import YTFieldNotFound


@requires_module("firefly")
def test_firefly_JSON_string():

    ds = fake_particle_ds()
    ad = ds.all_data()
    reader = ad.create_firefly_object(
        None,
        velocity_units="cm/s",
        coordinate_units="cm",
    )

    reader.writeToDisk(write_to_disk=False, file_extension=".json")

    ## reader.JSON was not output to string correctly
    ##  either Firefly is damaged or needs a hotfix-- try reinstalling.
    ##  if that doesn't work contact the developers
    ##  at github.com/ageller/Firefly/issues.
    assert len(reader.JSON) > 0


@requires_module("firefly")
def test_firefly_write_to_disk():
    tmpdir = tempfile.mkdtemp()

    ds = fake_particle_ds()
    ad = ds.all_data()
    reader = ad.create_firefly_object(
        tmpdir,
        velocity_units="cm/s",
        coordinate_units="cm",
    )

    reader.writeToDisk()


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
    ds_negative = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    ds = fake_particle_ds(
        fields=ds_fields,
        units=ds_field_units,
        negative=ds_negative,
    )
    return ds


@requires_module("firefly")
@pytest.mark.parametrize(
    "fields_to_include,fields_units,pgs_to_test",
    [
        (None, None, None),  # Test default values
        ([], [], None),  # Test empty fields
        (
            ["Masses"],
            ["code_mass"],
            ["PartType1", "gas"],
        ),  # Test common field (Masses)
        (
            ["Temperature"],
            ["code_temperature"],
            "gas",
        ),  # Test unique field (Temperature)
    ],
)
def test_field_string_specification(
    firefly_test_dataset, fields_to_include, fields_units, pgs_to_test
):
    dd = firefly_test_dataset.all_data()
    reader = dd.create_firefly_object(
        fields_to_include=fields_to_include,
        fields_units=fields_units,
        coordinate_units="code_length",
    )
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
        for idx, ptype in enumerate(all_pgs):
            pg = reader.particleGroups[idx]
            if ptype not in pgs_to_test:
                assert field not in pg.field_names
            else:
                assert field in pg.field_names
                arrind = np.flatnonzero(pg.field_names == field)[0]
                assert_array_equal(dd[ptype, field].d, pg.field_arrays[arrind])


@requires_module("firefly")
@pytest.mark.parametrize(
    "fields_to_include,fields_units",
    [
        (
            [("gas", "Temperature")],
            ["code_temperature"],
        ),  # Test existing field tuple (gas, Temperature)
        (
            [("PartType1", "Masses")],
            ["code_mass"],
        ),  # Test that tuples only bring in referenced particleGroup
    ],
)
def test_field_tuple_specification(
    firefly_test_dataset,
    fields_to_include,
    fields_units,
):
    dd = firefly_test_dataset.all_data()
    reader = dd.create_firefly_object(
        fields_to_include=fields_to_include,
        fields_units=fields_units,
        coordinate_units="code_length",
    )
    assert_array_equal(
        dd[("PartType1", "relative_particle_position")].d,
        reader.particleGroups[0].coordinates,
    )
    assert_array_equal(
        dd[("gas", "relative_particle_position")].d,
        reader.particleGroups[1].coordinates,
    )
    all_pgs = reader.particleGroups
    all_pgs_names = ["PartType1", "gas"]
    for field in fields_to_include:
        ftype, fname = field
        for pgi in range(2):
            pg = all_pgs[pgi]
            if ftype == all_pgs_names[pgi]:
                assert fname in pg.field_names
                arrind = np.flatnonzero(pg.field_names == fname)[0]
                assert_array_equal(dd[field].d, pg.field_arrays[arrind])
            else:
                assert fname not in pg.field_names


@requires_module("firefly")
@pytest.mark.parametrize(
    "fields_to_include,fields_units",
    [
        (
            ["dinos"],
            ["code_length"],
        ),  # Test nonexistent field (dinos)
        (
            [("PartType1", "Temperature")],
            ["code_temperature"],
        ),  # Test nonexistent field tuple (PartType1, Temperature)
    ],
)
def test_invalid_field_specifications(
    firefly_test_dataset, fields_to_include, fields_units
):

    dd = firefly_test_dataset.all_data()
    with pytest.raises(YTFieldNotFound):
        dd.create_firefly_object(
            fields_to_include=fields_to_include,
            fields_units=fields_units,
            coordinate_units="code_length",
        )


@requires_module("firefly")
def test_field_mixed_specification(firefly_test_dataset):
    dd = firefly_test_dataset.all_data()

    reader = dd.create_firefly_object(
        fields_to_include=["Masses", ("gas", "Temperature")],
        fields_units=["code_mass", "code_temperature"],
    )

    PartType1 = reader.particleGroups[0]
    gas = reader.particleGroups[1]
    assert "Masses" in PartType1.field_names
    arrind = np.flatnonzero(PartType1.field_names == "Masses")[0]
    assert_array_equal(dd[("PartType1", "Masses")].d, PartType1.field_arrays[arrind])
    assert "Masses" in gas.field_names
    arrind = np.flatnonzero(gas.field_names == "Masses")[0]
    assert_array_equal(dd[("gas", "Masses")].d, gas.field_arrays[arrind])

    assert "Temperature" not in PartType1.field_names
    assert "Temperature" in gas.field_names
    arrind = np.flatnonzero(gas.field_names == "Temperature")[0]
    assert_array_equal(dd[("gas", "Temperature")].d, gas.field_arrays[arrind])
