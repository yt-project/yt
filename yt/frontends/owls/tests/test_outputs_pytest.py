"""
Title: test_owls.py
Purpose: OWLS frontend tests using the snapshot_033 dataset
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.data_objects.particle_filters import add_particle_filter
from yt.frontends.owls.api import OWLSDataset
from yt.testing import ParticleSelectionComparison
from yt.utilities.answer_testing.answer_tests import sph_answer

# Test data
os33 = "snapshot_033/snap_033.0.hdf5"


@pytest.mark.answer_test
class TestOwls:
    answer_file = None

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [os33], indirect=True)
    def test_snapshot_033(self, f, w, d, a, ds):
        self.hashes.update(sph_answer(ds, "snap_033", 2 * 128 ** 3, f, w, d, a))

    @pytest.mark.big_data
    @pytest.mark.parametrize("ds", [os33], indirect=True)
    def test_owls_psc(self, ds):
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()

    @pytest.mark.parametrize("ds", [os33], indirect=True)
    def test_OWLSDataset(self, ds):
        assert isinstance(ds, OWLSDataset)

    @pytest.mark.parametrize("ds", [os33], indirect=True)
    def test_OWLS_particlefilter(self, ds):
        ad = ds.all_data()

        def cold_gas(pfilter, data):
            temperature = data[pfilter.filtered_type, "Temperature"]
            filter = temperature.in_units("K") <= 1e5
            return filter

        add_particle_filter(
            "gas_cold",
            function=cold_gas,
            filtered_type="PartType0",
            requires=["Temperature"],
        )
        ds.add_particle_filter("gas_cold")
        mask = ad["PartType0", "Temperature"] <= 1e5
        assert (
            ad["PartType0", "Temperature"][mask].shape
            == ad["gas_cold", "Temperature"].shape
        )
