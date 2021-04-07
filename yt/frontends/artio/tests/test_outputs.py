from yt.frontends.artio.api import ARTIODataset
from yt.loaders import load
from yt.testing import (
    assert_allclose_units,
    assert_equal,
    requires_file,
    units_override_check,
)
from yt.utilities.answer_testing.framework import (
    FieldValuesTest,
    PixelizedProjectionValuesTest,
    create_obj,
    data_dir_load,
    requires_ds,
)

_fields = (
    ("gas", "temperature"),
    ("gas", "density"),
    ("gas", "velocity_magnitude"),
    ("deposit", "all_density"),
    ("deposit", "all_count"),
)

sizmbhloz = "sizmbhloz-clref04SNth-rs9_a0.9011/sizmbhloz-clref04SNth-rs9_a0.9011.art"


@requires_ds(sizmbhloz)
def test_sizmbhloz():
    ds = data_dir_load(sizmbhloz)
    ds.max_range = 1024 * 1024
    assert_equal(str(ds), "sizmbhloz-clref04SNth-rs9_a0.9011.art")
    dso = [None, ("sphere", ("max", (0.1, "unitary")))]
    for dobj_name in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, ("gas", "density")]:
                    yield PixelizedProjectionValuesTest(
                        ds, axis, field, weight_field, dobj_name
                    )
            yield FieldValuesTest(ds, field, dobj_name)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj[("index", "ones")].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        assert_equal(s1, s2)
    assert_equal(ds.particle_type_counts, {"N-BODY": 100000, "STAR": 110650})


@requires_file(sizmbhloz)
def test_ARTIODataset():
    assert isinstance(data_dir_load(sizmbhloz), ARTIODataset)


@requires_file(sizmbhloz)
def test_units_override():
    units_override_check(sizmbhloz)


@requires_file(sizmbhloz)
def test_particle_derived_field():
    def star_age_alias(field, data):
        # test to make sure we get back data in the correct units
        # during field detection
        return data["STAR", "age"].in_units("Myr")

    ds = load(sizmbhloz)

    ds.add_field(
        ("STAR", "new_field"),
        function=star_age_alias,
        units="Myr",
        sampling_type="particle",
    )

    ad = ds.all_data()

    assert_allclose_units(ad["STAR", "age"].in_units("Myr"), ad["STAR", "new_field"])
