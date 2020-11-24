import pytest

from yt.frontends.gamer.api import GAMERDataset
from yt.testing import requires_file, units_override_check
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    grid_hierarchy,
    grid_values,
    parentage_relationships,
    projection_values,
)

# Test data
jet = "InteractingJets/jet_000002"
psiDM = "WaveDarkMatter/psiDM_000020"
plummer = "Plummer/plummer_000000"
mhd_vortex = "MHDOrszagTangVortex/Data_000018"

jet_units = {
    "length_unit": (1.0, "kpc"),
    "time_unit": (3.08567758096e13, "s"),
    "mass_unit": (1.4690033e36, "g"),
}

axes = [0, 1, 2]
objs = [None, ("sphere", ("max", (0.1, "unitary")))]
weights = [None, "density"]
jet_fields = ["temperature", "density", "velocity_magnitude"]
psi_fields = ["Dens", "Real", "Imag"]
plummer_fields = [("gamer", "ParDens"), ("deposit", "io_cic")]
mhd_fields = [("gamer", "CCMagX"), ("gamer", "CCMagY"), ("gas", "magnetic_energy")]
fields = [
    jet_fields,
    psi_fields,
    plummer_fields,
    mhd_fields,
]
ds_list = [
    [jet, {"kwargs": {"units_override": jet_units}}],
    psiDM,
    plummer,
    mhd_vortex,
]


def get_pairs():
    pairs = []
    for i, ds in enumerate(ds_list):
        for f in fields[i]:
            pairs.append((ds, f))
    return pairs


@pytest.mark.answer_test
class TestGamer:
    answer_file = None
    saved_hashes = None

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", ds_list, indirect=True)
    def test_grid_hierarchy_parentage_relationships(self, ds, big_data):
        if str(ds) == "rps.0062" and not big_data:
            pytest.skip("--answer-big-data not used.")
        self.hashes.update({"grid_hierarchy": grid_hierarchy(ds)})
        self.hashes.update({"parentage_relationships": parentage_relationships(ds)})

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    def test_grid_values(self, f, ds, big_data):
        if str(ds) == "rps.0062" and not big_data:
            pytest.skip("--answer-big-data not used.")
        self.hashes.update({"grid_values": grid_values(ds, f)})

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", objs, indirect=True)
    def test_field_values(self, d, f, ds, big_data):
        if str(ds) == "rps.0062" and not big_data:
            pytest.skip("--answer-big-data not used.")
        self.hashes.update({"field_values": field_values(ds, f, d)})

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", objs, indirect=True)
    @pytest.mark.parametrize("a", axes, indirect=True)
    @pytest.mark.parametrize("w", weights, indirect=True)
    def test_projection_values(self, a, d, w, f, ds, big_data):
        if str(ds) == "rps.0062" and not big_data:
            pytest.skip("--answer-big-data not used.")
        self.hashes.update({"projection_values": projection_values(ds, a, f, w, d)})

    @pytest.mark.parametrize("ds", [psiDM], indirect=True)
    def test_GAMERDataset(self, ds):
        assert isinstance(ds, GAMERDataset)

    @requires_file(jet)
    def test_units_override(self):
        units_override_check(jet)
