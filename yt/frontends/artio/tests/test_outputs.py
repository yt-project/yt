import pytest

from yt.frontends.artio.api import ARTIODataset
from yt.testing import assert_allclose_units, assert_equal, units_override_check
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    pixelized_projection_values,
)
from yt.utilities.answer_testing.testing_utilities import create_obj
from yt.utilities.answer_testing.testing_utilities import requires_ds

# Data file
sizmbhloz = "sizmbhloz-clref04SNth-rs9_a0.9011/"
sizmbhloz += "sizmbhloz-clref04SNth-rs9_a0.9011.art"

# Test data
axes = [0, 1, 2]
objs = [None, ("sphere", ("max", (0.1, "unitary")))]
weights = [None, "density"]

# Which velocity magnitude? gas? all? just velocity_magnitude isn't a field
# Does fv need particle_type?
fields = [
    "temperature",
    "density",
    "velocity_magnitude",
    ("deposit", "all_density"),
    ("deposit", "all_count"),
]


@pytest.mark.answer_test
class TestArtIo:
    answer_file = None
    saved_hashes = None

    @pytest.mark.parametrize("ds", [sizmbhloz], indirect=True)
    @pytest.mark.parametrize("d", objs, indirect=True)
    def test_sizmbhloz_validation(self, d, ds):
        ds.max_range = 1024 * 1024
        dobj = create_obj(ds, d)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        assert_equal(s1, s2)
        assert_equal(ds.particle_type_counts, {"N-BODY": 100000, "STAR": 110650})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [sizmbhloz], indirect=True)
    @pytest.mark.parametrize("d", objs, indirect=True)
    @pytest.mark.parametrize("a", axes, indirect=True)
    @pytest.mark.parametrize("w", weights, indirect=True)
    @pytest.mark.parametrize("f", fields, indirect=True)
    def test_sizmbhloz_pixelized_projection_values(self, d, a, w, f, ds):
        ds.max_range = 1024 * 1024
        ppv = pixelized_projection_values(ds, a, f, w, d)
        self.hashes.update({"pixelized_projection_values": ppv})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [sizmbhloz], indirect=True)
    @pytest.mark.parametrize("d", objs, indirect=True)
    @pytest.mark.parametrize("f", fields, indirect=True)
    def test_sizmbhloz_field_values(self, d, f, ds):
        ds.max_range = 1024 * 1024
        fv = field_values(ds, f, d)
        self.hashes.update({"field_values": fv})

    @pytest.mark.parametrize("ds", [sizmbhloz], indirect=True)
    def test_ARTIODataset(self, ds):
        assert isinstance(ds, ARTIODataset)

    @requires_ds(sizmbhloz)
    def test_units_override(self):
        units_override_check(sizmbhloz)

    @pytest.mark.parametrize("ds", [sizmbhloz], indirect=True)
    def test_particle_derived_field(self, ds):
        def star_age_alias(field, data):
            # test to make sure we get back data in the correct units
            # during field detection
            return data["STAR", "age"].in_units("Myr")

        ds.add_field(
            ("STAR", "new_field"),
            function=star_age_alias,
            units="Myr",
            sampling_type="particle",
        )
        ad = ds.all_data()
        assert_allclose_units(
            ad["STAR", "age"].in_units("Myr"), ad["STAR", "new_field"]
        )
