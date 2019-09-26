"""
Title: test_gadget_fof.py
Purpose: GadgetFOF frontend tests using gadget_fof datasets
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from collections import OrderedDict

import numpy as np
import pytest

from yt.frontends.gadget_fof.api import GadgetFOFDataset
from yt.testing import \
    assert_array_equal, \
    assert_equal, \
    requires_file
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

# Test data
g5 = "gadget_fof_halos/groups_005/fof_subhalo_tab_005.0.hdf5"
g42 = "gadget_fof_halos/groups_042/fof_subhalo_tab_042.0.hdf5"
g298 = "gadget_halos/data/groups_298/fof_subhalo_tab_298.0.hdf5"
g56 = "gadget_halos/data/groups_056/fof_subhalo_tab_056.0.hdf5"
g76 = "gadget_halos/data/groups_076/fof_subhalo_tab_076.0.hdf5"


#============================================
#               TestGadgetFOF
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestGadgetFOF(fw.AnswerTest):
    #-----
    # test_fields_g5
    #-----
    @requires_file(g5)
    def test_fields_g5(self, ds_g5):
        fields = ("particle_position_x", "particle_position_y",
                   "particle_position_z", "particle_velocity_x",
                   "particle_velocity_y", "particle_velocity_z",
                   "particle_mass", "particle_identifier")
        hashes = OrderedDict()
        hashes['field_values'] = OrderedDict()
        for field in fields:
            fv_hd = utils.generate_hash(
                self.field_values_test(ds_g5, field, particle_type=True)
            )
            hashes['field_values'][field] = fv_hd
        hashes = {'fields_g5' : hashes}
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)

    #-----
    # test_fields_g42
    #-----
    @utils.requires_ds(g42)
    def test_fields_g42(self, ds_g42):
        fields = ("particle_position_x", "particle_position_y",
                   "particle_position_z", "particle_velocity_x",
                   "particle_velocity_y", "particle_velocity_z",
                   "particle_mass", "particle_identifier")
        hashes = OrderedDict()
        hashes['field_values'] = OrderedDict()
        for field in fields:
            fv_hd = utils.generate_hash(
                self.field_values_test(ds_g42, field, particle_type=True)
            )
            hashes['field_values'][field] = fv_hd
        hashes = {'fields_g42' : hashes}
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)

    #-----
    # test_GadgetFOFDataset
    #-----
    @requires_file(g42)
    def test_GadgetFOFDataset(self, ds_g42):
        assert isinstance(ds_g42, GadgetFOFDataset)

    #-----
    # test_subhalos
    #-----
    @requires_file(g298)
    def test_subhalos(self, ds_g298):
        total_sub = 0
        total_int = 0
        for hid in range(0, ds_g298.index.particle_count["Group"]):
            my_h = ds_g298.halo("Group", hid)
            h_ids = my_h["ID"]
            for sid in range(int(my_h["subhalo_number"][0])):
                my_s = ds_g298.halo("Subhalo", (my_h.particle_identifier, sid))
                total_sub += my_s["ID"].size
                total_int += np.intersect1d(h_ids, my_s["ID"]).size
        assert_equal(total_sub, total_int)

    #-----
    # test_halo_masses
    #-----
    @requires_file(g298)
    def test_halo_masses(self, ds_g298):
        ad = ds_g298.all_data()
        for ptype in ["Group", "Subhalo"]:
            nhalos = ds_g298.index.particle_count[ptype]
            mass = ds_g298.arr(np.zeros(nhalos), "code_mass")
            for i in range(nhalos):
                halo = ds_g298.halo(ptype, i)
                mass[i] = halo.mass
            assert_array_equal(ad[ptype, "particle_mass"], mass)

    #-----
    # test_unbalanced_dataset
    #-----
    @requires_file(g56)
    def test_unbalanced_dataset(self, ds_g56):
        halo = ds_g56.halo("Group", 0)
        halo["member_ids"]
        assert True

    #-----
    # test_3file_halo
    #-----
    @requires_file(g76)
    def test_3file_halo(self, ds_g76):
        # this halo's particles are distributed over 3 files with the
        # middle file being empty
        halo = ds_g76.halo("Group", 6)
        halo["member_ids"]
        assert True
