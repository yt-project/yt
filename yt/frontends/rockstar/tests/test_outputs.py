import pytest

from yt.frontends.rockstar.api import RockstarDataset
from yt.testing import ParticleSelectionComparison
from yt.utilities.answer_testing.answer_tests import field_values

# Test data
r1 = "rockstar_halos/halos_0.0.bin"


@pytest.mark.answer_test
class TestRockstar:
    answer_file = None
    saved_hashes = None

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [r1], indirect=True)
    def test_fields_r1(self, f, ds):
        fv = field_values(ds, f, particle_type=True)
        self.hashes.update({"field_values": fv})

    @pytest.mark.parametrize("ds", [r1], indirect=True)
    def test_RockstarDataset(self, ds):
        assert isinstance(ds, RockstarDataset)

    @pytest.mark.parametrize("ds", [r1], indirect=True)
    def test_particle_selection(self, ds):
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()
