import pytest

from yt.testing import fake_random_ds
from yt.utilities.answer_testing.answer_tests import extract_connected_sets


@pytest.mark.answer_test
class TestConnectedSets:
    self.answer_file = None

    @pytest.mark.usefixtures("hashing")
    def test_connected_sets(self):
        ds = fake_random_ds(16, nprocs=8, particles=16 ** 3)
        data_source = ds.disk([0.5, 0.5, 0.5], [0.0, 0.0, 1.0], (8, "kpc"), (1, "kpc"))
        field = ("gas", "density")
        min_val, max_val = data_source[field].min() / 2, data_source[field].max() / 2
        data_source.extract_connected_sets(
            field, 5, min_val, max_val, log_space=True, cumulative=True
        )
        ecs = extract_connected_sets(ds, data_source, field, 5, min_val, max_val)
        self.hashes.update({"extract_connected_sets": ecs})
