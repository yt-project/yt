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
from yt.utilities.answer_testing import utils

# Test data
ahf_halos = 'ahf_halos/snap_N64L16_135.parameter'


@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestAHF:
    @requires_file(ahf_halos)
    def test_AHFHalosDataset(self, ds_ahf_halos):
        assert isinstance(ds_ahf_halos, AHFHalosDataset)
        ad = ds_ahf_halos.all_data()
        ad['particle_mass']
        psc = ParticleSelectionComparison(ds_ahf_halos)
        psc.run_defaults()

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(ahf_halos)
    def test_fields_ahf_halos(self, field, ds_ahf_halos):
        fv = field_values(ds_ahf_halos, field, particle_type=True)
        self.hashes.update({'field_values' : fv})
