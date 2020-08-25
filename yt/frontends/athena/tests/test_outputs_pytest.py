import pytest

import yt.units as u
from yt.frontends.athena.api import AthenaDataset
from yt.loaders import load
from yt.testing import (
    assert_allclose_units,
    assert_almost_equal,
    assert_equal,
    disable_dataset_cache,
    requires_file,
)
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    grid_hierarchy,
    grid_values,
    parentage_relationships,
    projection_values,
)

# Test data
cloud = "ShockCloud/id0/Cloud.0050.vtk"
blast = "MHDBlast/id0/Blast.0100.vtk"
stripping = "RamPressureStripping/id0/rps.0062.vtk"
sloshing = "MHDSloshing/virgo_low_res.0054.vtk"


# Test data params
uo_stripping = {
    "time_unit": 3.086e14,
    "length_unit": 8.0236e22,
    "mass_unit": 9.999e-30 * 8.0236e22 ** 3,
}
uo_blast = {
    "length_unit": (1.0, "pc"),
    "mass_unit": (2.38858753789e-24, "g/cm**3*pc**3"),
    "time_unit": (1.0, "s*pc/km"),
}
uo_sloshing = {
    "length_unit": (1.0, "Mpc"),
    "time_unit": (1.0, "Myr"),
    "mass_unit": (1.0e14, "Msun"),
}


ds_list = [
    cloud,
    blast,
    [stripping, {"kwargs": {"units_override": uo_stripping}}],
]
a_list = [0, 1, 2]
w_list = [None, "density"]
d_list = [None, ("sphere", ("max", (0.1, "unitary")))]
cloud_fields = [
    ("athena", "scalar[0]"),
    ("athena", "density"),
    ("athena", "total_energy"),
]
blast_fields = [
    ("gas", "temperature"),
    ("athena", "density"),
    ("gas", "velocity_magnitude"),
]
stripping_fields = [
    ("gas", "temperature"),
    ("athena", "density"),
    ("athena", "specific_scalar[0]"),
]
f_list = [
    cloud_fields,
    blast_fields,
    stripping_fields,
]


def get_pairs():
    pairs = []
    for i, ds in enumerate(ds_list):
        for f in f_list[i]:
            pairs.append((ds, f))
    return pairs


@pytest.mark.answer_test
class TestAthena:
    answer_file = None

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", ds_list, indirect=True)
    def test_gh_pr(self, ds, big_data):
        if str(ds) == "rps.0062" and not big_data:
            pytest.skip("--answer-big-data not used.")
        self.hashes.update({"grid_hierarchy": grid_hierarchy(ds)})
        self.hashes.update({"parentage_relationships": parentage_relationships(ds)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    def test_gv(self, f, ds, big_data):
        if str(ds) == "rps.0062" and not big_data:
            pytest.skip("--answer-big-data not used.")
        self.hashes.update({"grid_values": grid_values(ds, f)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", d_list, indirect=True)
    def test_fv(self, d, f, ds, big_data):
        if str(ds) == "rps.0062" and not big_data:
            pytest.skip("--answer-big-data not used.")
        self.hashes.update({"field_values": field_values(ds, f, d)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", d_list, indirect=True)
    @pytest.mark.parametrize("a", a_list, indirect=True)
    @pytest.mark.parametrize("w", w_list, indirect=True)
    def test_pv(self, a, d, w, f, ds, big_data):
        if str(ds) == "rps.0062" and not big_data:
            pytest.skip("--answer-big-data not used.")
        self.hashes.update({"projection_values": projection_values(ds, a, f, w, d)})

    @requires_file(blast)
    def test_blast_override(self):
        # verify that overriding units causes derived unit values to be updated.
        # see issue #1259
        ds = load(blast, units_override=uo_blast)
        assert_almost_equal(
            float(ds.magnetic_unit.in_units("gauss")), 5.478674679698131e-07, decimal=14
        )

    @requires_file(sloshing)
    @disable_dataset_cache
    def test_nprocs(self):
        ds1 = load(sloshing, units_override=uo_sloshing)
        sp1 = ds1.sphere("c", (100.0, "kpc"))
        prj1 = ds1.proj("density", 0)
        ds2 = load(sloshing, units_override=uo_sloshing, nprocs=8)
        sp2 = ds2.sphere("c", (100.0, "kpc"))
        prj2 = ds1.proj("density", 0)
        ds3 = load(sloshing, parameters=uo_sloshing)
        assert_equal(ds3.length_unit, u.Mpc)
        assert_equal(ds3.time_unit, u.Myr)
        assert_equal(ds3.mass_unit, 1e14 * u.Msun)
        assert_equal(
            sp1.quantities.extrema("pressure"), sp2.quantities.extrema("pressure")
        )
        assert_allclose_units(
            sp1.quantities.total_quantity("pressure"),
            sp2.quantities.total_quantity("pressure"),
        )
        for ax in "xyz":
            assert_equal(
                sp1.quantities.extrema(f"velocity_{ax}"),
                sp2.quantities.extrema(f"velocity_{ax}"),
            )
        assert_allclose_units(
            sp1.quantities.bulk_velocity(), sp2.quantities.bulk_velocity()
        )
        assert_equal(prj1["density"], prj2["density"])

    @pytest.mark.parametrize("ds", [cloud], indirect=True)
    def test_AthenaDataset(self, ds):
        assert isinstance(ds, AthenaDataset)
