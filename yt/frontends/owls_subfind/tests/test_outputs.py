import pytest

from yt.utilities.answer_testing.answer_tests import field_values

# Test data
g1 = "owls_fof_halos/groups_001/group_001.0.hdf5"
g8 = "owls_fof_halos/groups_008/group_008.0.hdf5"

g_fields = [
    "particle_position_x",
    "particle_position_y",
    "particle_position_z",
    "particle_mass",
]


@pytest.mark.answer_test
class TestOwlsSubfind:
    answer_file = None
    saved_hashes = None

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [g1, g8], indirect=True)
    @pytest.mark.parametrize("f", g_fields, indirect=True)
    def test_g_fields(self, f, ds):
        fv = field_values(ds, f, particle_type=True)
        self.hashes.update({"field_values": fv})
