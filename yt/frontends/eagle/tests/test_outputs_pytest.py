import pytest

from yt.frontends.eagle.api import EagleDataset
from yt.testing import ParticleSelectionComparison

# Test data
s28 = "snapshot_028_z000p000/snap_028_z000p000.0.hdf5"
s399 = "snipshot_399_z000p000/snip_399_z000p000.0.hdf5"


@pytest.mark.answer_test
class TestEagle:
    answer_file = None
    saved_hashes = None

    @pytest.mark.parametrize("ds", [s28], indirect=True)
    def test_EagleDataset(self, ds):
        assert isinstance(ds, EagleDataset)
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()

    @pytest.mark.parametrize("ds", [s399], indirect=True)
    def test_Snipshot(self, ds):
        assert isinstance(ds, EagleDataset)
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()
