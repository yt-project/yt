"""
Title: test_ahf.py
Purpose: Contains functions for answer testing the AHF frontend
Notes:
    Copyright (c) 2017, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.frontends.ahf.api import AHFHalosDataset
from yt.testing import ParticleSelectionComparison, requires_file
from yt.utilities.answer_testing.answer_tests import field_values

# Test data
ahf_halos = "ahf_halos/snap_N64L16_135.parameter"
ahf_kwargs = {"kwargs" : {"hubble_constant": 0.7}}

ahf_fields = [
    ("nbody", "particle_position_x"),
    ("nbody", "particle_position_y"),
    ("nbody", "particle_position_z"),
    ("nbody", "particle_mass"),
]


@pytest.mark.answer_test
class TestAHF:
    @requires_file(ahf_halos)
    @pytest.mark.parametrize("ds", [[ahf_halos, ahf_kwargs]], indirect=True)
    def test_AHFHalosDataset(self, ds):
        assert isinstance(ds, AHFHalosDataset)
        ad = ds.all_data()
        ad["particle_mass"]
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [[ahf_halos, ahf_kwargs]], indirect=True)
    @pytest.mark.parametrize("f", ahf_fields, indirect=True)
    def test_fields_ahf_halos(self, f, ds):
        fv = field_values(ds, f, particle_type=True)
        self.hashes.update({"field_values": fv})
