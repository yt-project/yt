"""
Title: test_gadget_fof.py
Purpose: GadgetFOF frontend tests using gadget_fof datasets
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import numpy as np
import pytest

from yt.frontends.gadget_fof.api import GadgetFOFDataset
from yt.testing import \
    assert_array_equal, \
    assert_equal, \
    requires_file
    ParticleSelectionComparison
from yt.utilities.answer_testing.answer_tests import field_values
from yt.utilities.answer_testing.utils import requires_ds, data_dir_load

# Test data
g5 = "gadget_fof_halos/groups_005/fof_subhalo_tab_005.0.hdf5"
g42 = "gadget_fof_halos/groups_042/fof_subhalo_tab_042.0.hdf5"
g298 = "gadget_halos/data/groups_298/fof_subhalo_tab_298.0.hdf5"
g56 = "gadget_halos/data/groups_056/fof_subhalo_tab_056.0.hdf5"
g76 = "gadget_halos/data/groups_076/fof_subhalo_tab_076.0.hdf5"


@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestGadgetFOF:

    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [g5], indirect=True)
    def test_fields_g5(self, field, ds):
        fv = field_values(ds, field, particle_type=True)
        self.hashes.update({'field_values' : fv})

    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [g42], indirect=True)
    def test_fields_g42(self, field, ds):
        fv = field_values(ds, field, particle_type=True)
        self.hashes.update({'field_values' : fv})

    @pytest.mark.parametrize('ds', [g42], indirect=True)
    def test_GadgetFOFDataset(self, ds):
        assert isinstance(ds, GadgetFOFDataset)

    @pytest.mark.parametrize('ds', [g298], indirect=True)
    def test_subhalos(self, ds):
        total_sub = 0
        total_int = 0
        for hid in range(0, ds.index.particle_count["Group"]):
            my_h = ds.halo("Group", hid)
            h_ids = my_h["ID"]
            for sid in range(int(my_h["subhalo_number"][0])):
                my_s = ds.halo("Subhalo", (my_h.particle_identifier, sid))
                total_sub += my_s["ID"].size
                total_int += np.intersect1d(h_ids, my_s["ID"]).size
        assert_equal(total_sub, total_int)

    @pytest.mark.parametrize('ds', [g298], indirect=True)
    def test_halo_masses(self, ds):
        ad = ds.all_data()
        for ptype in ["Group", "Subhalo"]:
            nhalos = ds.index.particle_count[ptype]
            mass = ds.arr(np.zeros(nhalos), "code_mass")
            for i in range(nhalos):
                halo = ds.halo(ptype, i)
                mass[i] = halo.mass
            assert_array_equal(ad[ptype, "particle_mass"], mass)

    @pytest.mark.parametrize('ds', [g56], indirect=True)
    def test_unbalanced_dataset(self, ds):
        halo = ds.halo("Group", 0)
        assert_equal(len(halo["member_ids"]), 33)
        assert_equal(halo["member_ids"].min().d, 723254.0)
        assert_equal(halo["member_ids"].max().d, 772662.0)

    @pytest.mark.parametrize('ds', [g76], indirect=True)
    def test_3file_halo(self, ds):
        # this halo's particles are distributed over 3 files with the
        # middle file being empty
        halo = ds.halo("Group", 6)
        halo["member_ids"]
        assert True

    @pytest.mark.parametrize('ds', [g298], indirect=True)
    def test_particle_selection(self, ds):
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()
