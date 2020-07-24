"""
Title: test_owls_subfind.py
Purpose: OWLSSubfind frontend tests using owls_fof_halos datasets
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.utilities.answer_testing.answer_tests import field_values

# Test data
g1 = "owls_fof_halos/groups_001/group_001.0.hdf5"
g8 = "owls_fof_halos/groups_008/group_008.0.hdf5"


@pytest.mark.answer_test
@pytest.mark.usefixtures("answer_file")
class TestOwlsSubfind:
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [g8], indirect=True)
    def test_fields_g8(self, field, ds):
        fv = field_values(ds, field, particle_type=True)
        self.hashes.update({"field_values": fv})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [g1], indirect=True)
    def test_fields_g1(self, field, ds):
        fv = field_values(ds, field, particle_type=True)
        self.hashes.update({"field_values": fv})
