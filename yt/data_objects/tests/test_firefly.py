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
        ("pt1", "particle_position_x"),
        ("pt1", "particle_position_y"),
        ("pt1", "particle_position_z"),
        ("pt2", "particle_position_x"),
        ("pt2", "particle_position_y"),
        ("pt2", "particle_position_z"),
        # User input
        ("pt1", "common_field"),
        ("pt2", "common_field"),
        ("pt2", "pt2only_field"),
    ]
    ds_field_units = ["code_length"] * 9
    ds_negative = [0] * 9
    ds = fake_particle_ds(
        fields=ds_fields,
        units=ds_field_units,
        negative=ds_negative,
    )
    return ds


@requires_module("firefly")
@pytest.mark.parametrize(
    "fields_to_include,fields_units",
    [
        (None, None),  # Test default values
        ([], []),  # Test empty fields
    ],
)
def test_field_empty_specification(
    firefly_test_dataset, fields_to_include, fields_units
):
    dd = firefly_test_dataset.all_data()
    reader = dd.create_firefly_object(
        fields_to_include=fields_to_include,
        fields_units=fields_units,
        coordinate_units="code_length",
    )
    assert_array_equal(
        dd[("pt1", "relative_particle_position")].d,
        reader.particleGroups[0].coordinates,
    )
    assert_array_equal(
        dd[("pt2", "relative_particle_position")].d,
        reader.particleGroups[1].coordinates,
    )


@requires_module("firefly")
def test_field_string_specification(firefly_test_dataset):
    # Test unique field (pt2only_field)
    dd = firefly_test_dataset.all_data()
    # Unique field string will fallback to "all" field type ...
    with pytest.raises(YTFieldNotFound):
        reader = dd.create_firefly_object(
            fields_to_include=["pt2only_field"],
            fields_units=["code_length"],
            coordinate_units="code_length",
        )
        pytest.fail()

    # ... unless field type has been previously explicitly
    # specified
    dd["pt2", "pt2only_field"]
    reader = dd.create_firefly_object(
        fields_to_include=["pt2only_field"],
        fields_units=["code_length"],
        coordinate_units="code_length",
    )

    pt1 = reader.particleGroups[0]
    pt2 = reader.particleGroups[1]
    assert_array_equal(
        dd[("pt1", "relative_particle_position")].d,
        pt1.coordinates,
    )
    assert_array_equal(
        dd[("pt2", "relative_particle_position")].d,
        pt2.coordinates,
    )
    assert "pt2only_field" not in pt1.field_names
    assert "pt2only_field" in pt2.field_names
    arrind = np.flatnonzero(pt2.field_names == "pt2only_field")[0]
    assert_array_equal(dd[("pt2", "pt2only_field")].d, pt2.field_arrays[arrind])


@requires_module("firefly")
@pytest.mark.parametrize(
    "fields_to_include,fields_units",
    [
        (
            [("pt2", "pt2only_field")],
            ["code_length"],
        ),  # Test existing field tuple (pt2, pt2only_field)
        (
            [("pt1", "common_field")],
            ["code_length"],
        ),  # Test that tuples only bring in referenced particleGroup
        (
            [("all", "common_field")],
            ["code_length"],
        ),  # Test that "all" brings in all particleGroups
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
        dd[("pt1", "relative_particle_position")].d,
        reader.particleGroups[0].coordinates,
    )
    assert_array_equal(
        dd[("pt2", "relative_particle_position")].d,
        reader.particleGroups[1].coordinates,
    )
    all_pgs = reader.particleGroups
    all_pgs_names = ["pt1", "pt2"]
    for field in fields_to_include:
        ftype, fname = field
        for pgi in range(2):
            pg = all_pgs[pgi]
            if ftype == all_pgs_names[pgi]:
                assert fname in pg.field_names
                arrind = np.flatnonzero(pg.field_names == fname)[0]
                assert_array_equal(dd[field].d, pg.field_arrays[arrind])
            elif ftype == "all":
                assert fname in pg.field_names
                this_pg_name = all_pgs_names[pgi]
                arrind = np.flatnonzero(pg.field_names == fname)[0]
                assert_array_equal(dd[this_pg_name, fname].d, pg.field_arrays[arrind])
            else:
                assert fname not in pg.field_names


@requires_module("firefly")
@pytest.mark.parametrize(
    "fields_to_include,fields_units,ErrorType",
    [
        (
            ["dinos"],
            ["code_length"],
            YTFieldNotFound,
        ),  # Test nonexistent field (dinos)
        (
            ["common_field"],
            ["code_length"],
            ValueError,
        ),  # Test ambiguous field (common_field)
        (
            [("pt1", "pt2only_field")],
            ["code_length"],
            YTFieldNotFound,
        ),  # Test nonexistent field tuple (pt1, pt2only_field)
    ],
)
def test_field_invalid_specification(
    firefly_test_dataset, fields_to_include, fields_units, ErrorType
):

    dd = firefly_test_dataset.all_data()
    with pytest.raises(ErrorType):
        dd.create_firefly_object(
            fields_to_include=fields_to_include,
            fields_units=fields_units,
            coordinate_units="code_length",
        )


@requires_module("firefly")
def test_field_mixed_specification(firefly_test_dataset):
    dd = firefly_test_dataset.all_data()

    # Assume particle type has already been specified
    dd["pt2", "pt2only_field"]
    reader = dd.create_firefly_object(
        fields_to_include=["pt2only_field", ("pt1", "common_field")],
        fields_units=["code_length", "code_length"],
    )

    pt1 = reader.particleGroups[0]
    pt2 = reader.particleGroups[1]
    assert "common_field" in pt1.field_names
    assert "common_field" not in pt2.field_names
    arrind = np.flatnonzero(pt1.field_names == "common_field")[0]
    assert_array_equal(dd[("pt1", "common_field")].d, pt1.field_arrays[arrind])

    assert "pt2only_field" not in pt1.field_names
    assert "pt2only_field" in pt2.field_names
    arrind = np.flatnonzero(pt2.field_names == "pt2only_field")[0]
    assert_array_equal(dd[("pt2", "pt2only_field")].d, pt2.field_arrays[arrind])
